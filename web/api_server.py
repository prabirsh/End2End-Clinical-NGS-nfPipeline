#!/usr/bin/env python3
"""
Clinical NGS Pipeline - Web API Server
Maintained by: Prabir

A Flask-based API server for managing NGS analysis runs through a web interface.
Provides endpoints for job submission, status monitoring, and result retrieval.
"""

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import subprocess
import os
import json
from datetime import datetime
from pathlib import Path
import threading
import time

app = Flask(__name__)
CORS(app)

# Configuration
BASE_DIR = Path(__file__).parent.parent
WORK_DIR = BASE_DIR / "work"
OUTPUT_DIR = BASE_DIR / "results"
WORK_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

# Job tracking
jobs = {}
job_counter = 0

class PipelineJob:
    """
    Represents a pipeline analysis job with status tracking
    Maintained by: Prabir
    """
    def __init__(self, job_id, params):
        self.job_id = job_id
        self.params = params
        self.status = "queued"
        self.start_time = None
        self.end_time = None
        self.progress = 0
        self.logs = []
        self.results = {}
        
    def add_log(self, level, message):
        """Add a log entry"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.logs.append({
            "timestamp": timestamp,
            "level": level,
            "message": message
        })
        
    def to_dict(self):
        """Convert job to dictionary"""
        return {
            "job_id": self.job_id,
            "status": self.status,
            "params": self.params,
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "end_time": self.end_time.isoformat() if self.end_time else None,
            "progress": self.progress,
            "logs": self.logs[-50:],  # Last 50 logs
            "results": self.results
        }

def run_pipeline(job):
    """
    Execute the Nextflow pipeline
    Maintained by: Prabir
    """
    try:
        job.status = "running"
        job.start_time = datetime.now()
        job.add_log("info", f"Starting {job.params['analysis_type']} analysis")
        job.add_log("info", "Variant Caller: GATK HaplotypeCaller (Germline) / Mutect2 (Somatic)")
        
        # Construct Nextflow command
        cmd = [
            "nextflow", "run", str(BASE_DIR / "main.nf"),
            "-profile", job.params.get("profile", "docker"),
            "--analysis_type", job.params["analysis_type"],
            "--sample_sheet", job.params["sample_sheet"],
            "--reference_genome", job.params["reference_genome"],
            "--target_bed", job.params["target_bed"],
            "--outdir", str(OUTPUT_DIR / f"job_{job.job_id}"),
            "-work-dir", str(WORK_DIR / f"job_{job.job_id}")
        ]
        
        # Add Panel of Normals if tumor-only
        if job.params["analysis_type"] == "somatic_tumor_only" and "pon" in job.params:
            cmd.extend(["--pon", job.params["pon"]])
        
        job.add_log("info", f"Command: {' '.join(cmd)}")
        
        # Simulate pipeline execution (replace with actual Nextflow call in production)
        stages = [
            ("FastQC", "Quality control completed", 10),
            ("BWA-MEM", "Alignment completed", 30),
            ("Picard", "Duplicate marking completed", 20),
            ("GATK", "Base quality score recalibration completed", 15),
            ("Variant Calling", f"Variant calling completed using GATK {'HaplotypeCaller' if 'germline' in job.params['analysis_type'] else 'Mutect2'}", 50),
            ("VEP", "Annotation completed", 75),
            ("Report", "Clinical report generated", 100)
        ]
        
        for stage_name, message, progress in stages:
            time.sleep(2)  # Simulate processing time
            job.progress = progress
            job.add_log("success", f"{stage_name}: {message}")
        
        # In production, replace above simulation with:
        # process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # for line in process.stdout:
        #     job.add_log("info", line.strip())
        # process.wait()
        
        job.status = "completed"
        job.end_time = datetime.now()
        job.progress = 100
        job.add_log("success", "Pipeline completed successfully!")
        
        # Collect results
        job.results = {
            "vcf_file": f"job_{job.job_id}/variants/filtered.vcf.gz",
            "bam_file": f"job_{job.job_id}/alignment/aligned.bam",
            "report": f"job_{job.job_id}/reports/clinical_report.html",
            "qc_report": f"job_{job.job_id}/qc/multiqc_report.html"
        }
        
    except Exception as e:
        job.status = "failed"
        job.end_time = datetime.now()
        job.add_log("error", f"Pipeline failed: {str(e)}")

@app.route('/')
def index():
    """Serve the web interface"""
    return send_from_directory(BASE_DIR / "web", "index.html")

@app.route('/api/jobs', methods=['POST'])
def create_job():
    """Create a new pipeline job"""
    global job_counter
    
    data = request.json
    required_fields = ['analysis_type', 'sample_sheet', 'reference_genome', 'target_bed']
    
    # Validate required fields
    for field in required_fields:
        if field not in data:
            return jsonify({"error": f"Missing required field: {field}"}), 400
    
    # Validate analysis type
    valid_types = ['germline', 'somatic_paired', 'somatic_tumor_only']
    if data['analysis_type'] not in valid_types:
        return jsonify({"error": f"Invalid analysis_type. Must be one of: {', '.join(valid_types)}"}), 400
    
    # Check for PoN requirement
    if data['analysis_type'] == 'somatic_tumor_only' and 'pon' not in data:
        return jsonify({"error": "Panel of Normals (pon) is required for tumor-only analysis"}), 400
    
    # Create job
    job_counter += 1
    job = PipelineJob(job_counter, data)
    jobs[job_counter] = job
    
    # Start pipeline in background thread
    thread = threading.Thread(target=run_pipeline, args=(job,))
    thread.daemon = True
    thread.start()
    
    return jsonify({"job_id": job_counter, "status": "queued"}), 201

@app.route('/api/jobs/<int:job_id>', methods=['GET'])
def get_job(job_id):
    """Get job status and details"""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404
    
    return jsonify(jobs[job_id].to_dict())

@app.route('/api/jobs', methods=['GET'])
def list_jobs():
    """List all jobs"""
    return jsonify([job.to_dict() for job in jobs.values()])

@app.route('/api/jobs/<int:job_id>/logs', methods=['GET'])
def get_job_logs(job_id):
    """Get job logs"""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404
    
    return jsonify({"logs": jobs[job_id].logs})

@app.route('/api/jobs/<int:job_id>/results', methods=['GET'])
def get_job_results(job_id):
    """Get job results"""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404
    
    job = jobs[job_id]
    if job.status != "completed":
        return jsonify({"error": "Job not completed yet"}), 400
    
    return jsonify(job.results)

@app.route('/api/jobs/<int:job_id>/cancel', methods=['POST'])
def cancel_job(job_id):
    """Cancel a running job"""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404
    
    job = jobs[job_id]
    if job.status in ["completed", "failed", "cancelled"]:
        return jsonify({"error": f"Job already {job.status}"}), 400
    
    job.status = "cancelled"
    job.end_time = datetime.now()
    job.add_log("warning", "Job cancelled by user")
    
    return jsonify({"message": "Job cancelled successfully"})

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        "status": "healthy",
        "version": "2.0",
        "maintained_by": "Prabir",
        "timestamp": datetime.now().isoformat()
    })

@app.route('/api/info', methods=['GET'])
def pipeline_info():
    """Get pipeline information"""
    return jsonify({
        "name": "Clinical NGS Pipeline",
        "version": "2.0",
        "maintained_by": "Prabir",
        "variant_callers": {
            "germline": "GATK HaplotypeCaller 4.5+",
            "somatic": "GATK Mutect2 4.5+ with Panel of Normals"
        },
        "supported_analyses": [
            "germline",
            "somatic_paired",
            "somatic_tumor_only"
        ],
        "features": [
            "GATK HaplotypeCaller for germline variants",
            "GATK Mutect2 for somatic variants",
            "Panel of Normals support",
            "Comprehensive QC and coverage analysis",
            "Clinical-grade variant annotation",
            "Professional HTML reports"
        ]
    })

if __name__ == '__main__':
    print("="*60)
    print("Clinical NGS Pipeline - Web API Server")
    print("Maintained by: Prabir")
    print("="*60)
    print(f"Base Directory: {BASE_DIR}")
    print(f"Output Directory: {OUTPUT_DIR}")
    print(f"Work Directory: {WORK_DIR}")
    print("="*60)
    print("Starting server on http://0.0.0.0:5000")
    print("="*60)
    
    app.run(host='0.0.0.0', port=5000, debug=True)
