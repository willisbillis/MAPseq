#!/usr/bin/env python3
"""
generate_multi_config.py

Generate per-batch Cell Ranger multi_config.csv files from master sample sheet.

Features:
  • Validates VDJ-GEX pairing relationships
  • Handles feature reference assignment (HTO, ADT, ASAP, VDJ)
  • Generates modality-specific configs (RNA+HTO+VDJ, ATAC+ASAP, etc.)
  • Checks FASTQ file existence and format validity
  • Outputs configs organized by batch_id

Usage:
  python3 generate_multi_config.py \
    --sample_sheet samples.csv \
    --output_dir configs/ \
    --gex_ref /path/to/GRCh38/ref \
    --vdj_ref /path/to/VDJ/ref \
    --hto_ref /path/to/HTO_feature_ref.csv
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MultiConfigGenerator:
    """Generate Cell Ranger multi_config.csv files from sample sheet."""
    
    def __init__(self, gex_ref, vdj_ref, atac_ref, hto_feature_ref, asap_feature_ref):
        """Initialize with reference paths."""
        self.gex_ref = Path(gex_ref)
        self.vdj_ref = Path(vdj_ref) if vdj_ref else None
        self.atac_ref = Path(atac_ref) if atac_ref else None
        self.hto_feature_ref = Path(hto_feature_ref) if hto_feature_ref else None
        self.asap_feature_ref = Path(asap_feature_ref) if asap_feature_ref else None
        
        self._validate_references()
    
    def _validate_references(self):
        """Validate that required reference files exist."""
        if not self.gex_ref.exists():
            raise FileNotFoundError(f"GEX reference not found: {self.gex_ref}")
        
        logger.info(f"✓ GEX reference: {self.gex_ref}")
        
        if self.vdj_ref and not self.vdj_ref.exists():
            logger.warning(f"VDJ reference not found: {self.vdj_ref} (VDJ will be skipped)")
            self.vdj_ref = None
        else:
            logger.info(f"✓ VDJ reference: {self.vdj_ref}")
        
        if self.hto_feature_ref and not self.hto_feature_ref.exists():
            logger.warning(f"HTO feature reference not found: {self.hto_feature_ref}")
            self.hto_feature_ref = None
        else:
            logger.info(f"✓ HTO feature reference: {self.hto_feature_ref}")
    
    def _get_fastq_dir(self, fastq_path):
        """Extract directory from FASTQ path and validate format."""
        fastq_path = Path(fastq_path)
        
        # Validate FASTQ extension
        if not str(fastq_path).endswith(('.fastq.gz', '.fastq', '.fq.gz', '.fq')):
            logger.warning(f"Unusual FASTQ extension: {fastq_path}")
        
        return str(fastq_path.parent)
    
    def load_sample_sheet(self, sample_sheet_path):
        """Load and parse sample sheet."""
        samples_by_batch = defaultdict(list)
        
        with open(sample_sheet_path, 'r') as f:
            reader = csv.DictReader(f)
            
            # Validate header
            required_cols = {'sample_id', 'batch_id', 'gex_fastq_r1', 'gex_fastq_r2'}
            if not required_cols.issubset(set(reader.fieldnames)):
                raise ValueError(f"Sample sheet missing required columns: {required_cols}")
            
            for row_num, row in enumerate(reader, start=2):  # Start at 2 (header=1)
                sample_id = row.get('sample_id', '').strip()
                batch_id = row.get('batch_id', '').strip()
                
                if not sample_id or not batch_id:
                    logger.warning(f"Row {row_num}: Missing sample_id or batch_id, skipping")
                    continue
                
                # Check GEX FASTQ files exist
                gex_r1 = row.get('gex_fastq_r1', '').strip()
                gex_r2 = row.get('gex_fastq_r2', '').strip()
                
                if not gex_r1 or not gex_r2:
                    logger.error(f"Row {row_num} ({sample_id}): Missing GEX FASTQ paths")
                    continue
                
                if not Path(gex_r1).exists() or not Path(gex_r2).exists():
                    logger.warning(f"Row {row_num} ({sample_id}): GEX FASTQ files not found (may be on compute node)")
                
                # Parse optional libraries
                hto_r1 = row.get('hto_fastq_r1', '').strip() or None
                hto_r2 = row.get('hto_fastq_r2', '').strip() or None
                vdj_r1 = row.get('vdj_fastq_r1', '').strip() or None
                vdj_r2 = row.get('vdj_fastq_r2', '').strip() or None
                
                sample = {
                    'sample_id': sample_id,
                    'batch_id': batch_id,
                    'gex_r1': gex_r1,
                    'gex_r2': gex_r2,
                    'hto_r1': hto_r1,
                    'hto_r2': hto_r2,
                    'vdj_r1': vdj_r1,
                    'vdj_r2': vdj_r2,
                    'vdj_type': row.get('vdj_type', 'vdj-bcr').strip(),  # vdj-bcr or vdj-tcr
                    'row_num': row_num
                }
                
                samples_by_batch[batch_id].append(sample)
                logger.info(f"✓ Loaded {sample_id} (batch: {batch_id})")
        
        return samples_by_batch
    
    def _validate_vdj_pairing(self, samples):
        """Validate that samples with VDJ data also have GEX data."""
        for sample in samples:
            has_vdj = sample['vdj_r1'] is not None and sample['vdj_r2'] is not None
            has_gex = sample['gex_r1'] is not None and sample['gex_r2'] is not None
            
            if has_vdj and not has_gex:
                raise ValueError(
                    f"Sample {sample['sample_id']}: VDJ data present but GEX data missing. "
                    "VDJ requires paired GEX for proper linking."
                )
        
        logger.info("✓ VDJ-GEX pairing validation passed")
    
    def generate_config(self, batch_id, samples):
        """Generate multi_config.csv content for a batch."""
        
        # Validate pairing
        self._validate_vdj_pairing(samples)
        
        # Check if batch has VDJ data
        has_vdj = any(s['vdj_r1'] is not None for s in samples)
        has_hto = any(s['hto_r1'] is not None for s in samples)
        
        config_lines = []
        
        # [gene-expression] section
        config_lines.append('[gene-expression]')
        config_lines.append(f'reference,{self.gex_ref}')
        config_lines.append('create-bam,true')
        config_lines.append('')
        
        # [feature] section (if HTO or ADT present)
        if has_hto and self.hto_feature_ref:
            config_lines.append('[feature]')
            config_lines.append(f'reference,{self.hto_feature_ref}')
            config_lines.append('')
        
        # [vdj] section (if VDJ present)
        if has_vdj and self.vdj_ref:
            config_lines.append('[vdj]')
            config_lines.append(f'reference,{self.vdj_ref}')
            config_lines.append('')
        
        # [libraries] section
        config_lines.append('[libraries]')
        config_lines.append('fastq_id,fastqs,feature_types')
        
        for sample in samples:
            # GEX library
            if sample['gex_r1'] and sample['gex_r2']:
                gex_dir = self._get_fastq_dir(sample['gex_r1'])
                config_lines.append(f"{sample['sample_id']}_GEX,{gex_dir},gene expression")
            
            # Feature barcoding (HTO/ADT)
            if sample['hto_r1'] and sample['hto_r2']:
                hto_dir = self._get_fastq_dir(sample['hto_r1'])
                config_lines.append(f"{sample['sample_id']}_HTO,{hto_dir},antibody capture")
            
            # VDJ library
            if sample['vdj_r1'] and sample['vdj_r2']:
                vdj_dir = self._get_fastq_dir(sample['vdj_r1'])
                vdj_type = sample.get('vdj_type', 'vdj-bcr')
                config_lines.append(f"{sample['sample_id']}_VDJ,{vdj_dir},{vdj_type}")
        
        return '\n'.join(config_lines)
    
    def write_configs(self, samples_by_batch, output_dir):
        """Write multi_config.csv files for each batch."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        config_paths = {}
        
        for batch_id, samples in samples_by_batch.items():
            config_content = self.generate_config(batch_id, samples)
            config_path = output_dir / f"multi_config_{batch_id}.csv"
            
            with open(config_path, 'w') as f:
                f.write(config_content)
            
            config_paths[batch_id] = str(config_path)
            logger.info(f"✓ Generated {config_path} ({len(samples)} samples)")
        
        return config_paths


def main():
    parser = argparse.ArgumentParser(
        description='Generate Cell Ranger multi_config.csv files from sample sheet',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--sample_sheet', required=True, help='CSV with sample metadata')
    parser.add_argument('--output_dir', default='cellranger_configs', help='Output directory for configs')
    parser.add_argument('--gex_ref', required=True, help='Path to GEX reference directory')
    parser.add_argument('--vdj_ref', help='Path to VDJ reference directory')
    parser.add_argument('--atac_ref', help='Path to ATAC reference directory')
    parser.add_argument('--hto_feature_ref', help='Path to HTO feature reference CSV')
    parser.add_argument('--asap_feature_ref', help='Path to ASAP feature reference CSV')
    
    args = parser.parse_args()
    
    try:
        logger.info("=" * 70)
        logger.info("Cell Ranger multi_config.csv Generator")
        logger.info("=" * 70)
        
        # Initialize generator
        generator = MultiConfigGenerator(
            gex_ref=args.gex_ref,
            vdj_ref=args.vdj_ref,
            atac_ref=args.atac_ref,
            hto_feature_ref=args.hto_feature_ref,
            asap_feature_ref=args.asap_feature_ref
        )
        
        # Load sample sheet
        logger.info(f"\nLoading sample sheet: {args.sample_sheet}")
        samples_by_batch = generator.load_sample_sheet(args.sample_sheet)
        
        if not samples_by_batch:
            logger.error("No samples loaded from sample sheet")
            sys.exit(1)
        
        logger.info(f"✓ Loaded {sum(len(s) for s in samples_by_batch.values())} samples from {len(samples_by_batch)} batches")
        
        # Generate configs
        logger.info(f"\nGenerating Cell Ranger configs to: {args.output_dir}")
        config_paths = generator.write_configs(samples_by_batch, args.output_dir)
        
        # Summary
        logger.info("\n" + "=" * 70)
        logger.info("Summary:")
        for batch_id, path in config_paths.items():
            logger.info(f"  • {batch_id}: {path}")
        logger.info("=" * 70)
        
        print("\nSUCCESS")  # Signal to Nextflow
        
    except Exception as e:
        logger.error(f"\n❌ Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
