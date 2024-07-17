# How to view QC reports

### 1: Run the workflow with personalized report generation configurations

### 2: Navigate to qc directory:

```bash
cd ./qc
```
### 3: Create a Python web server:

```bash
python -m http.server
```
### 4: Open the link generated in a browser and open qc_report.html


### Note:
 Can view the whole slice images and the flatfield corrections without running the web server. However, to view the volume rendered brain, the web server is required.