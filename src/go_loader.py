"""
GO (Gene Ontology) annotation loader.
Parses GOA GAF files and maps to proteins.
"""

import os
import gzip
import logging
from typing import Dict, Set, List, Optional
from collections import defaultdict

logger = logging.getLogger(__name__)


class GOLoader:
    """
    Load GO annotations from GOA GAF files.
    """
    
    def __init__(self, cache_dir: str = "cache"):
        """
        Initialize GO loader.
        
        Args:
            cache_dir: Directory for caching GAF files
        """
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
    
    def load_from_gaf(self, gaf_file: str, taxid: Optional[int] = None, 
                     use_symbol: bool = False) -> Dict[str, Set[str]]:
        """
        Load GO annotations from GAF file.
        
        Args:
            gaf_file: Path to GAF file (can be .gz)
            taxid: Optional taxonomy ID to filter (if None, loads all)
            use_symbol: If True, use DB_Object_Symbol (column 2) instead of DB_Object_ID (column 1)
                       This is needed for SGD GAF files where symbols are yeast ORF names
            
        Returns:
            protein_go_terms: Dict mapping protein ID to set of GO term IDs
        """
        protein_go_terms = defaultdict(set)
        
        logger.info(f"Loading GO annotations from {gaf_file}...")
        if use_symbol:
            logger.info("Using DB_Object_Symbol for protein IDs (SGD format)")
        
        # Try different encodings
        encodings = ['utf-8', 'latin-1', 'iso-8859-1', 'cp1252']
        
        for encoding in encodings:
            try:
                if gaf_file.endswith('.gz'):
                    open_func = lambda f, e=encoding: gzip.open(f, 'rt', encoding=e)
                else:
                    open_func = lambda f, e=encoding: open(f, 'r', encoding=e, errors='replace')
                
                with open_func(gaf_file) as f:
                    # Test read
                    test_line = f.readline()
                    if test_line:
                        # Reset and process
                        f.seek(0)
                        break
            except (UnicodeDecodeError, AttributeError, TypeError):
                continue
        else:
            # Fallback: use errors='replace'
            if gaf_file.endswith('.gz'):
                open_func = lambda f: gzip.open(f, 'rt', encoding='utf-8', errors='replace')
            else:
                open_func = lambda f: open(f, 'r', encoding='utf-8', errors='replace')
        
        with open_func(gaf_file) as f:
            for line in f:
                if line.startswith('!'):
                    continue  # Skip comment lines
                
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                
                # GAF format fields
                # Column 0: DB
                # Column 1: DB_Object_ID
                # Column 2: DB_Object_Symbol (e.g., YDL159W for SGD)
                # Column 3: Qualifier
                # Column 4: GO_ID
                
                if use_symbol:
                    protein_id = parts[2] if len(parts) > 2 else parts[1]  # Use symbol (SGD format)
                else:
                    protein_id = parts[1]  # Use DB_Object_ID (GOA format)
                
                qualifier = parts[3] if len(parts) > 3 else ""
                go_id = parts[4] if len(parts) > 4 else None
                
                if not go_id or not go_id.startswith('GO:'):
                    continue
                
                # Get taxon_id if available (usually column 12 in GAF 2.1)
                taxon_id = None
                if len(parts) > 12:
                    taxon_id = parts[12]
                
                # Filter by taxid if provided
                if taxid is not None:
                    if taxon_id:
                        taxon_parts = taxon_id.split('|')
                        if not any(str(taxid) in tp for tp in taxon_parts):
                            continue
                    else:
                        # If taxid filtering requested but no taxon_id found, skip
                        continue
                
                # Only include 'NOT' excluded terms and valid evidence codes
                if qualifier == 'NOT':
                    continue
                
                # Filter by evidence code (optional - include all for now)
                # Common codes: EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, IBA, IBD, IKR, IRD, ISS, ISO, ISA, ISM, IGC, RCA, TAS, IC, ND
                
                if protein_id and go_id:
                    protein_go_terms[protein_id].add(go_id)
        
        logger.info(f"Loaded GO annotations for {len(protein_go_terms)} proteins")
        return dict(protein_go_terms)
    
    def get_go_terms_for_cluster(self, cluster_proteins: Set[str], 
                                  protein_go_terms: Dict[str, Set[str]]) -> Set[str]:
        """
        Get all GO terms associated with proteins in a cluster.
        
        Args:
            cluster_proteins: Set of protein IDs in cluster
            protein_go_terms: Dict mapping protein ID to GO terms
            
        Returns:
            Set of GO term IDs
        """
        cluster_go_terms = set()
        for protein in cluster_proteins:
            if protein in protein_go_terms:
                cluster_go_terms.update(protein_go_terms[protein])
        return cluster_go_terms

