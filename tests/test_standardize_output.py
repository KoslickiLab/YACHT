import unittest 
import os
import subprocess

class TestScript(unittest.TestCase):
    def test_everything_exists(self):
        script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        script_dir = os.path.join(script_dir, 'srcs')
        script_full_path = os.path.join(script_dir, 'standardize_yacht_output.py')
        assert os.path.exists(script_full_path)

        yacht_output = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/result.xlsx')
        assert os.path.exists(yacht_output)

        genome_to_taxid = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid.tsv')
        assert os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata')
        assert os.path.exists(outdir)

        cmd = f"python {script_full_path} --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert os.path.exists(os.path.join(outdir, 'cami_result.cami'))

    def test_wrong_yacht_output(self):
        script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        script_dir = os.path.join(script_dir, 'srcs')
        script_full_path = os.path.join(script_dir, 'standardize_yacht_output.py')
        assert os.path.exists(script_full_path)

        yacht_output = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/result_nonexisting.xlsx')
        assert not os.path.exists(yacht_output)

        genome_to_taxid = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid.tsv')
        assert os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata')
        assert os.path.exists(outdir)
        
        cmd = f"python {script_full_path} --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        with self.assertRaises(subprocess.CalledProcessError):
            res = subprocess.run(cmd, shell=True, check=True)

    def test_wrong_genome_to_taxid(self):
        script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        script_dir = os.path.join(script_dir, 'srcs')
        script_full_path = os.path.join(script_dir, 'standardize_yacht_output.py')
        assert os.path.exists(script_full_path)

        yacht_output = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/result.xlsx')
        assert os.path.exists(yacht_output)

        genome_to_taxid = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid_nonexisting.tsv')
        assert not os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata')
        assert os.path.exists(outdir)

        cmd = f"python {script_full_path} --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        with self.assertRaises(subprocess.CalledProcessError):
            res = subprocess.run(cmd, shell=True, check=True)

    def test_wrong_outdir(self):
        script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        script_dir = os.path.join(script_dir, 'srcs')
        script_full_path = os.path.join(script_dir, 'standardize_yacht_output.py')
        assert os.path.exists(script_full_path)

        yacht_output = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/result.xlsx')
        assert os.path.exists(yacht_output)

        genome_to_taxid = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid.tsv')
        assert os.path.exists(genome_to_taxid)

        outdir = os.path.join(os.path.dirname(__file__), 'testdata_nonexisting')
        cmd = 'rm -rf ' + outdir
        try:
            subprocess.run(cmd, shell=True, check=True)
        except:
            pass
        assert not os.path.exists(outdir)

        cmd = f"python {script_full_path} --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert os.path.exists(outdir)

        cmd = 'rm -rf ' + outdir
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0


if __name__ == '__main__':
    unittest.main()