import numpy as np
import os
from Bio import SeqIO
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import f1_score
from .taxonomy import TaxonomyCache, TaxonomyLabeler

def read_fasta_folder(folder_path):
    sequences = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith(('.fasta', '.fa', '.fna')):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r') as handle:
                first_record = next(SeqIO.parse(handle, "fasta"), None)
                if first_record:
                    file_id = os.path.splitext(filename)[0]
                    sequences.append((file_id, str(first_record.seq)))
    
    return sequences

def lzset(sequence):
    dictionary = set()
    i = 0
    while i < len(sequence):
        j = 1
        while i + j <= len(sequence) and sequence[i:i+j] in dictionary:
            j += 1
        if i + j <= len(sequence):
            dictionary.add(sequence[i:i+j])
        i += max(1, j-1)
    return dictionary

def lzjd_distance(seq1, seq2):
    set1 = lzset(seq1)
    set2 = lzset(seq2)
    
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    return 1 - (intersection / union if union > 0 else 0)

def calculate_distance_matrix(sequences):
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    filenames = [seq[0] for seq in sequences]  

    for i in range(n):
        for j in range(i+1, n):
            dist = lzjd_distance(sequences[i][1], sequences[j][1])
            distance_matrix[i][j] = dist
            distance_matrix[j][i] = dist
            
    return distance_matrix, filenames

def generate_tree(distance_matrix, labels, method='average'):
    condensed_dist = squareform(distance_matrix)

    linkage_matrix = hierarchy.linkage(condensed_dist, method=method)
    
    def to_newick(node, labels, linkage_matrix):
        if node < len(labels):
            return labels[node]
        else:
            left = int(linkage_matrix[node-len(labels), 0])
            right = int(linkage_matrix[node-len(labels), 1])
            dist = linkage_matrix[node-len(labels), 2] / 2
            return f"({to_newick(left, labels, linkage_matrix)}:{dist:.4f},{to_newick(right, labels, linkage_matrix)}:{dist:.4f})"

    newick = to_newick(2*len(labels)-2, labels, linkage_matrix) + ";"
    return newick

def save_results(distance_matrix, filenames, output_prefix, tree_method='average'):
    matrix_file = f"{output_prefix}_matrix.txt"
    with open(matrix_file, 'w') as f:
        f.write('\t' + '\t'.join(filenames) + '\n')
        for i, filename in enumerate(filenames):
            row = [filename] + [str(round(dist, 4)) for dist in distance_matrix[i]]
            f.write('\t'.join(row) + '\n')

    tree_file = f"{output_prefix}_tree.nwk"
    try:
        newick = generate_tree(distance_matrix, filenames, method=tree_method)
        with open(tree_file, 'w') as f:
            f.write(newick)
    except Exception as e:
        print(f"Error: {e}")

class LZJDKNNCorrector:
    def __init__(self, sequences, taxonomy_level='genus', confidence_threshold=0.6):
        self.sequences = sequences
        self.sequence_names = [seq[0] for seq in sequences]
        self.taxonomy_level = taxonomy_level
        self.confidence_threshold = confidence_threshold
        self.taxonomy_labeler = TaxonomyLabeler()
        
        self.labels = [self.taxonomy_labeler.get_label(name, taxonomy_level) 
                      for name in self.sequence_names]

        self.distance_matrix = self._calculate_distance_matrix()
        
    def _calculate_distance_matrix(self):
        n = len(self.sequences)
        matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                dist = lzjd_distance(self.sequences[i][1], self.sequences[j][1])
                matrix[i,j] = matrix[j,i] = dist
                
        return matrix
    
    def perform_knn_classification(self, k=3):
        X = self.distance_matrix
        y = np.array(self.labels)

        unique_labels = sorted(list(set(self.labels)))
        n_classes = len(unique_labels)
        
        loo = LeaveOneOut()
        self.probability_matrix = np.zeros((len(X), n_classes))
        predictions = []
        true_labels = []
        
        label_to_idx = {label: idx for idx, label in enumerate(unique_labels)}
        
        knn = KNeighborsClassifier(n_neighbors=k, metric='precomputed')
        
        for train_idx, test_idx in loo.split(X):
            X_train = X[train_idx][:, train_idx]
            X_test = X[test_idx][:, train_idx]
            y_train = y[train_idx]
            
            knn.fit(X_train, y_train)
            
            pred_proba = knn.predict_proba(X_test)[0]
            
            full_prob = np.zeros(n_classes)
            for i, label in enumerate(knn.classes_):
                full_prob[label_to_idx[label]] = pred_proba[i]
            
            self.probability_matrix[test_idx] = full_prob
            
            predictions.append(knn.predict(X_test)[0])
            true_labels.append(y[test_idx][0])
        
        return predictions, true_labels
    
    def correct_distances(self, base_alpha=0.3):
        predictions, true_labels = self.perform_knn_classification()
        n = len(self.sequences)
        self.corrected_distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                prob_sim = np.dot(self.probability_matrix[i], self.probability_matrix[j])
                prob_sim /= (np.linalg.norm(self.probability_matrix[i]) * 
                           np.linalg.norm(self.probability_matrix[j]))
                
                tax_dist = self._calculate_taxonomic_distance(self.labels[i], self.labels[j])

                confidence = min(np.max(self.probability_matrix[i]),
                               np.max(self.probability_matrix[j]))
                
                is_correct = (predictions[i] == true_labels[i] and 
                            predictions[j] == true_labels[j])
                
                alpha = self._calculate_dynamic_alpha(base_alpha, confidence, is_correct)

                knn_dist = 0.5 * (1 - prob_sim) + 0.5 * tax_dist
                corrected_dist = (1 - alpha) * self.distance_matrix[i,j] + alpha * knn_dist
                
                self.corrected_distances[i,j] = self.corrected_distances[j,i] = corrected_dist
        
        return self.corrected_distances
    
    def _calculate_taxonomic_distance(self, label1, label2):
        if label1 == label2:
            return 0.0
            
        tax1 = self.taxonomy_labeler.get_taxonomy_info(label1)
        tax2 = self.taxonomy_labeler.get_taxonomy_info(label2)
        
        if not tax1 or not tax2:
            return 1.0
        
        hierarchy = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
        total_levels = len(hierarchy)
        
        for level in hierarchy:
            if (level in tax1 and level in tax2 and 
                tax1[level] == tax2[level]):
                return hierarchy.index(level) / total_levels
                
        return 1.0
    
    def _calculate_dynamic_alpha(self, base_alpha, confidence, is_correct):
        if not is_correct:
            confidence *= 0.5
        if confidence < self.confidence_threshold:
            return 0
        return base_alpha * confidence

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
            
            linkage_matrix = hierarchy.linkage(condensed_dist, method=tree_method)
            
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
