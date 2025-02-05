import os
import numpy as np
import re
from pathlib import Path
from multiprocessing import Pool
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import f1_score
from scipy.cluster.hierarchy import linkage
from Bio import Entrez

Entrez.email = "x1058513236@gmail.com"

from taxonomy import TaxonomyCache, TaxonomyLabeler

class BCDClassifier:
    def __init__(self, compressed_dir, taxonomy_level='class'):
        self.compressed_dir = Path(compressed_dir)
        self.compressed_sizes = {}
        self.distance_matrix = None
        self.labels = []
        self.sequence_names = []
        self.valid_extensions = ('.fasta', '.fa', '.fna')
        self.taxonomy_labeler = TaxonomyLabeler()
        self.taxonomy_level = taxonomy_level
        self.corrected_distances = None 
        self.knn_distances = None 
        self.confidence_threshold = 0.6
        
    def load_compressed_sizes(self):
        sequence_names = []
        labels = []

        for ext in self.valid_extensions:
            pattern = f"class_*{ext}.biox"
            for file_path in self.compressed_dir.glob(pattern):
                name = file_path.stem
                name = name.replace(ext, '')
                
                if name.startswith('class_'):
                    sequence_names.append(name)
                    label = self.taxonomy_labeler.get_label(name, self.taxonomy_level)
                    labels.append(label)
                    self.compressed_sizes[name] = file_path.stat().st_size

        for ext in self.valid_extensions:
            pattern = f"*-*{ext}.biox"
            for file_path in self.compressed_dir.glob(pattern):
                name = file_path.stem
                name = name.replace(ext, '')
                
                if not name.startswith('class_'):
                    self.compressed_sizes[name] = file_path.stat().st_size

        print(f"\nLoaded {len(sequence_names)} individual sequence files")
        print(f"Loaded {len(self.compressed_sizes) - len(sequence_names)} concatenated sequence files")

        self.sequence_names = sequence_names
        return sequence_names, labels

    def calculate_ncd(self, seq1_name, seq2_name):
        c_x = self.compressed_sizes[seq1_name]
        c_y = self.compressed_sizes[seq2_name]
            
        name1 = seq1_name.replace('class_', '')
        name2 = seq2_name.replace('class_', '')
            
        concat_name1 = f"{name1}-{name2}"
        concat_name2 = f"{name2}-{name1}"
            
        c_xy = self.compressed_sizes.get(concat_name1)
            
        if c_xy is None:
            c_xy = self.compressed_sizes.get(concat_name2)
            if c_xy is None:
                raise ValueError(f"Concatenated sequence file not found: {concat_name1} or {concat_name2}")
        
        ncd1 = (c_xy - c_x) / c_y
        ncd2 = (c_xy - c_y) / c_x
        
        w1 = c_y / (c_x + c_y) 
        w2 = c_x / (c_x + c_y)
        
        ncd = w1 * ncd1 + w2 * ncd2
        
        return ncd

    
    def build_distance_matrix(self, sequence_names, n_jobs=-1):
        n = len(sequence_names)
        self.distance_matrix = np.zeros((n, n))
        
        tasks = []
        for i in range(n):
            for j in range(i+1, n):
                tasks.append((sequence_names[i], sequence_names[j]))
        
        with Pool(processes=n_jobs) as pool:
            results = pool.starmap(self.calculate_ncd, tasks)
        
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                self.distance_matrix[i,j] = results[idx]
                self.distance_matrix[j,i] = results[idx]
                idx += 1
        
    
    def perform_knn_classification(self, k=3):
        X = self.distance_matrix
        y = np.array(self.labels)
        
        knn = KNeighborsClassifier(n_neighbors=k, metric='precomputed')
        
        loo = LeaveOneOut()
        predictions = []
        true_labels = []
        probabilities = []
        
        all_classes = sorted(set(self.labels))
        n_classes = len(all_classes)
        
        for train_idx, test_idx in loo.split(X):
            X_train, X_test = X[train_idx][:, train_idx], X[test_idx][:, train_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            knn.fit(X_train, y_train)
            pred = knn.predict(X_test)
            
            prob = knn.predict_proba(X_test)[0]
            
            full_prob = np.zeros(n_classes)
            for i, class_label in enumerate(knn.classes_):
                idx = all_classes.index(class_label)
                full_prob[idx] = prob[i]
            
            predictions.append(pred[0])
            true_labels.append(y_test[0])
            probabilities.append(full_prob)
        
        macro_f1 = f1_score(true_labels, predictions, average='macro')
        
        print(f"Macro F1 Score: {macro_f1:.4f}")

        self.probability_matrix = np.array(probabilities)
        return predictions, true_labels, probabilities, macro_f1
    
    def save_results(self, predictions, true_labels, probabilities, macro_f1, output_dir, tree_method='single'):
        os.makedirs(output_dir, exist_ok=True)

        basic_info_file = os.path.join(output_dir, 'basic_info.txt')
        with open(basic_info_file, 'w') as f:
            f.write("NCBI Classification Analysis Report\n")
            f.write("=" * 80 + "\n\n")
            f.write("Basic Information\n")
            f.write("-" * 40 + "\n")
            f.write(f"Taxonomy Level: {self.taxonomy_level}\n")
            f.write(f"Total Samples: {len(true_labels)}\n")
            f.write(f"Number of Classes: {len(set(true_labels))}\n")
            
            total = len(true_labels)
            correct = sum(1 for t, p in zip(true_labels, predictions) if t == p)
            accuracy = correct / total * 100
            
            f.write(f"Classification Accuracy: {accuracy:.2f}%\n")
            f.write(f"Macro F1 Score: {macro_f1:.4f}\n")
        
        ncd_file = os.path.join(output_dir, 'ncd_values.txt')
        with open(ncd_file, 'w') as f:
            f.write("Sequence Pair NCD Values\n")
            f.write("-" * 40 + "\n")
            n = len(self.sequence_names)
            for i in range(n):
                for j in range(i+1, n):
                    name1 = self.sequence_names[i].replace('class_', '')
                    name2 = self.sequence_names[j].replace('class_', '')
                    ncd = self.distance_matrix[i,j]
                    f.write(f"{name1}-{name2:<50}{ncd:.4f}\n")
        
        dist_file = os.path.join(output_dir, 'knn_corrected_distances.txt')
        with open(dist_file, 'w') as f:
            f.write("Corrected Distance Matrix\n")
            f.write("-" * 40 + "\n")
            n = len(self.sequence_names)
            f.write("Sequence")
            for name in self.sequence_names:
                f.write(f"\t{name.replace('class_', '')}")
            f.write("\n")
            for i in range(n):
                f.write(f"{self.sequence_names[i].replace('class_', '')}")
                for j in range(n):
                    f.write(f"\t{self.corrected_distances[i,j]:.4f}")
                f.write("\n")

        tree_file = os.path.join(output_dir, 'phylogenetic_tree.nwk')
        try:
            condensed_dist = []
            n = len(self.sequence_names)
            for i in range(n):
                for j in range(i+1, n):
                    condensed_dist.append(self.corrected_distances[i,j])
            
            linkage_matrix = linkage(condensed_dist, method=tree_method)
            
            labels = [name.replace('class_', '') for name in self.sequence_names]
            
            def to_newick(node, labels, linkage_matrix):
                if node < len(labels):
                    return labels[node]
                else:
                    left = int(linkage_matrix[node-len(labels), 0])
                    right = int(linkage_matrix[node-len(labels), 1])
                    dist = linkage_matrix[node-len(labels), 2] / 2
                    return f"({to_newick(left, labels, linkage_matrix)}:{dist:.4f},{to_newick(right, labels, linkage_matrix)}:{dist:.4f})"
            
            with open(tree_file, 'w') as f:
                newick = to_newick(2*len(labels)-2, labels, linkage_matrix) + ";"
                f.write(newick)
                
        except Exception as e:
            print(f"Error generating phylogenetic tree: {e}")
        
        print(f"\nResults saved to:")
        print(f"1. Basic Information: {basic_info_file}")
        print(f"2. NCD Values: {ncd_file}")
        print(f"3. KNN Corrected Distances: {dist_file}")
        print(f"4. Phylogenetic Tree: {tree_file}")

    def calculate_knn_distance(self, prob1, prob2, label1, label2):
        prob_sim = np.dot(prob1, prob2) / (np.linalg.norm(prob1) * np.linalg.norm(prob2))
        
        if label1 == label2:
            tax_dist = 0.0
        else:
            tax1 = self.taxonomy_labeler.get_taxonomy_info(label1)
            tax2 = self.taxonomy_labeler.get_taxonomy_info(label2)
            
            if not tax1 or not tax2:
                return 1.0  
            
            hierarchy = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
            total_levels = len(hierarchy)
            
            common_level = total_levels
            for level in hierarchy:
                if level in tax1 and level in tax2 and tax1[level] == tax2[level]:
                    common_level = hierarchy.index(level)
                    break
            
            tax_dist = common_level / total_levels
        
        w1, w2 = 0.5, 0.5  
        knn_dist = w1 * (1 - prob_sim) + w2 * tax_dist
        return knn_dist
    
    def calculate_dynamic_alpha(self, base_alpha, confidence, is_correct):
        if not is_correct:
            confidence *= 0.5  
        if confidence < self.confidence_threshold:
            return 0 
        return base_alpha * confidence
    
    def correct_distances(self, base_alpha=0.3):
        n = len(self.sequence_names)
        self.corrected_distances = np.zeros((n, n))
        self.knn_distances = np.zeros((n, n))
        
        knn = KNeighborsClassifier(n_neighbors=1, metric='precomputed')
        predictions = []
        confidences = []
        
        for i in range(n):
            mask = np.ones(n, dtype=bool)
            mask[i] = False
            knn.fit(self.distance_matrix[mask][:, mask], 
                   np.array(self.labels)[mask])
            pred = knn.predict(self.distance_matrix[i:i+1, mask])
            prob = knn.predict_proba(self.distance_matrix[i:i+1, mask])
            predictions.append(pred[0])
            confidences.append(np.max(prob[0]))
        
        for i in range(n):
            for j in range(i+1, n):
                knn_dist = self.calculate_knn_distance(
                    self.probability_matrix[i],
                    self.probability_matrix[j],
                    self.labels[i],
                    self.labels[j]
                )
                
                is_correct_i = predictions[i] == self.labels[i]
                is_correct_j = predictions[j] == self.labels[j]
                
                confidence = min(confidences[i], confidences[j])
                
                dynamic_alpha = self.calculate_dynamic_alpha(
                    base_alpha, 
                    confidence,
                    is_correct_i and is_correct_j
                )
                
                corrected_dist = (1 - dynamic_alpha) * self.distance_matrix[i,j] + \
                               dynamic_alpha * knn_dist
                
                self.knn_distances[i,j] = self.knn_distances[j,i] = knn_dist
                self.corrected_distances[i,j] = self.corrected_distances[j,i] = corrected_dist
        
        return self.corrected_distances 