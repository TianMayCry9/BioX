from Bio import Entrez
import sqlite3
from datetime import datetime, timedelta
import re

Entrez.email = "x1058513236@gmail.com"

class TaxonomyCache:
    def __init__(self, cache_file='taxonomy_cache.db'):
        self.cache_file = cache_file
        self._init_db()

    def _init_db(self):
        with sqlite3.connect(self.cache_file) as conn:
            conn.execute('''CREATE TABLE IF NOT EXISTS taxonomy_cache
                          (species_name TEXT PRIMARY KEY,
                           kingdom TEXT,
                           phylum TEXT,
                           class TEXT,
                           order_name TEXT,
                           family TEXT,
                           genus TEXT,
                           species TEXT,
                           timestamp DATETIME)''')

    def get(self, species_name):
        with sqlite3.connect(self.cache_file) as conn:
            cursor = conn.execute('''SELECT * FROM taxonomy_cache 
                                   WHERE species_name = ?''', (species_name,))
            result = cursor.fetchone()
            if result:
                timestamp = datetime.strptime(result[-1], '%Y-%m-%d %H:%M:%S')
                if datetime.now() - timestamp < timedelta(days=30):
                    return dict(zip(['kingdom', 'phylum', 'class', 'order', 
                                   'family', 'genus', 'species'], result[1:-1]))
        return None

    def set(self, species_name, taxonomy_info):
        with sqlite3.connect(self.cache_file) as conn:
            conn.execute('''INSERT OR REPLACE INTO taxonomy_cache VALUES 
                          (?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                       (species_name,
                        taxonomy_info.get('kingdom'),
                        taxonomy_info.get('phylum'),
                        taxonomy_info.get('class'),
                        taxonomy_info.get('order'),
                        taxonomy_info.get('family'),
                        taxonomy_info.get('genus'),
                        taxonomy_info.get('species'),
                        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

class TaxonomyLabeler:
    def __init__(self):
        self.unknown_counter = 1
        self.cache = TaxonomyCache()

    def get_taxonomy_info(self, species_name):
        try:
            cached_info = self.cache.get(species_name)
            if cached_info:
                return cached_info

            handle = Entrez.esearch(db="taxonomy", term=species_name)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                return None

            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            if not records:
                return None

            taxonomy = records[0]
            lineage = {}
            
            for item in taxonomy.get('LineageEx', []):
                rank = item.get('Rank')
                if rank in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']:
                    lineage[rank] = item.get('ScientificName')
            
            lineage['species'] = taxonomy.get('ScientificName')
            
            self.cache.set(species_name, lineage)
            return lineage

        except Exception as e:
            print(f"Error getting taxonomy info for {species_name}: {e}")
            return None

    def get_label(self, species_name, taxonomy_level):
        species_name = re.sub(r'class_|\..*$', '', species_name)
        
        species_name = re.sub(r'(?<!^)(?=[A-Z])', ' ', species_name)
        
        tax_info = self.get_taxonomy_info(species_name)
        if not tax_info or taxonomy_level not in tax_info:
            label = f"Unknown{self.unknown_counter}"
            self.unknown_counter += 1
            return label

        return tax_info[taxonomy_level]