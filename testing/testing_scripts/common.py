"""
Common code for unit testing of rules generated with Snakemake 7.32.4.
"""

from pathlib import Path
import subprocess as sp
import os
import difflib

class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )  
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if str(f).startswith(".snakemake"):
                    continue
                if f in expected_files:
                    self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                elif str(f).startswith("logs/") or str(f).startswith("benchmarks") or str(f).startswith(".java/") or str(f).startswith("qc/resources") or str(f) == "qc/qc_report.html":
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_files(self, generated_file, expected_file):
        sp.check_output(["cmp", generated_file, expected_file])
        
        
class ImperfectOutputChecker(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        if(os.path.getsize(generated_file) and os.path.getsize(expected_file)):
            with open(generated_file, 'rb') as gen, open(expected_file, 'rb') as exp:
                total_similarity = []
                while True:
                    gen_content = gen.read(1024)
                    exp_content = exp.read(1024)
                    if not gen_content or not exp_content:
                        break
                    similarity_ratio = float(difflib.SequenceMatcher(None, gen_content, exp_content).ratio())
                    total_similarity.append(similarity_ratio)
            final_sim_score = sum(total_similarity)/len(total_similarity)
            print(final_sim_score)
            assert final_sim_score>=0.995, final_sim_score
        elif os.path.getsize(generated_file) != os.path.getsize(expected_file):
            raise ValueError("Files not equal")
        
