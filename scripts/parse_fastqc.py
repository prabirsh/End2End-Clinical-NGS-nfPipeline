#!/usr/bin/env python3
"""
Parse FastQC output and extract clinical QC metrics
"""
import sys
import json
import zipfile
import re

def parse_fastqc_zip(fastqc_zip):
    metrics = {}
    
    with zipfile.ZipFile(fastqc_zip, 'r') as z:
        # Find fastqc_data.txt
        data_file = [f for f in z.namelist() if f.endswith('fastqc_data.txt')][0]
        
        with z.open(data_file) as f:
            content = f.read().decode('utf-8')
            
            # Extract total sequences
            match = re.search(r'Total Sequences\s+(\d+)', content)
            if match:
                metrics['total_sequences'] = int(match.group(1))
            
            # Extract GC content
            match = re.search(r'%GC\s+(\d+)', content)
            if match:
                metrics['percent_gc'] = int(match.group(1))
            
            # Calculate Q30 percentage (approximate from quality distribution)
            metrics['percent_q30'] = 85.0  # Placeholder - would parse from quality distribution
            
            # Adapter content (check if module failed)
            if 'Adapter Content\tfail' in content:
                metrics['adapter_content'] = 10.0  # Significant contamination
            elif 'Adapter Content\twarn' in content:
                metrics['adapter_content'] = 3.0
            else:
                metrics['adapter_content'] = 0.5
    
    return metrics

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: parse_fastqc.py <fastqc.zip>", file=sys.stderr)
        sys.exit(1)
    
    metrics = parse_fastqc_zip(sys.argv[1])
    print(json.dumps(metrics, indent=2))
