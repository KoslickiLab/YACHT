import unittest
import os
import subprocess

PATH_TO_YACHT_OUTPUT = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/result.xlsx')
PATH_TO_GENOME_TO_TAXID = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid.tsv')

def assert_file_exists(file_path):
    assert os.path.exists(file_path)

def assert_file_not_exists(file_path):
    assert not os.path.exists(file_path)

def create_outdir(outdir):
    cmd = f'rm -rf {outdir}'
    try:
        subprocess.run(cmd, shell=True, check=True)
    except:
        pass
    assert not os.path.exists(outdir)

def cleanup_outdir(outdir):
    cmd = f'rm -rf {outdir}'
    res = subprocess.run(cmd, shell=True, check=True)
    assert res.returncode == 0

class TestScript(unittest.TestCase):

    def test_everything_exists(self):
        yacht_output = PATH_TO_YACHT_OUTPUT
        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        outdir = os.path.join(os.path.dirname(__file__), 'testdata')

        assert_file_exists(yacht_output)
        assert_file_exists(genome_to_taxid)
        assert_file_exists(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert_file_exists(os.path.join(outdir, 'cami_result.cami'))

        cmd = 'rm ' + os.path.join(outdir, 'cami_result.cami')
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0

    def test_biom(self):
        yacht_output = PATH_TO_YACHT_OUTPUT
        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        outdir = os.path.join(os.path.dirname(__file__), 'testdata')

        assert_file_exists(yacht_output)
        assert_file_exists(genome_to_taxid)
        assert_file_exists(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir} --mode biom"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert_file_exists(os.path.join(outdir, 'cami_result.biom'))

        # remove the generated file
        cmd = 'rm ' + os.path.join(outdir, 'cami_result.biom')
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0

    def test_graphplan(self):
        yacht_output = PATH_TO_YACHT_OUTPUT
        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        outdir = os.path.join(os.path.dirname(__file__), 'testdata')

        assert_file_exists(yacht_output)
        assert_file_exists(genome_to_taxid)
        assert_file_exists(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir} --mode graphplan"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert_file_exists(os.path.join(outdir, 'cami_result.graphplan.newick'))

        # remove the generated file
        cmd = 'rm ' + os.path.join(outdir, 'cami_result.graphplan.newick')
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0

    def test_all_mode(self):
        yacht_output = PATH_TO_YACHT_OUTPUT
        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        outdir = os.path.join(os.path.dirname(__file__), 'testdata')

        assert_file_exists(yacht_output)
        assert_file_exists(genome_to_taxid)
        assert_file_exists(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir} --mode all"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert_file_exists(os.path.join(outdir, 'cami_result.cami'))
        assert_file_exists(os.path.join(outdir, 'cami_result.biom'))
        assert_file_exists(os.path.join(outdir, 'cami_result.graphplan.newick'))

        # remove the generated files
        cmd = 'rm ' + os.path.join(outdir, 'cami_result.cami')
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0

        cmd = 'rm ' + os.path.join(outdir, 'cami_result.biom')
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0

        cmd = 'rm ' + os.path.join(outdir, 'cami_result.graphplan.newick')
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0

    def test_wrong_mode(self):
        yacht_output = PATH_TO_YACHT_OUTPUT
        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        outdir = os.path.join(os.path.dirname(__file__), 'testdata')

        assert_file_exists(yacht_output)
        assert_file_exists(genome_to_taxid)
        assert_file_exists(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir} --mode wrong_mode"
        with self.assertRaises(subprocess.CalledProcessError):
            _ = subprocess.run(cmd, shell=True, check=True)

    def test_wrong_yacht_output(self):
        yacht_output = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/result_nonexisting.xlsx')
        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        outdir = os.path.join(os.path.dirname(__file__), 'testdata')

        assert_file_not_exists(yacht_output)
        assert_file_exists(genome_to_taxid)
        assert_file_exists(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        with self.assertRaises(subprocess.CalledProcessError):
            _ = subprocess.run(cmd, shell=True, check=True)

    def test_wrong_genome_to_taxid(self):
        yacht_output = PATH_TO_YACHT_OUTPUT
        genome_to_taxid = os.path.join(os.path.dirname(__file__), 'testdata/standardize_output_testdata/toy_genome_to_taxid_nonexisting.tsv')
        outdir = os.path.join(os.path.dirname(__file__), 'testdata')

        assert_file_exists(yacht_output)
        assert_file_not_exists(genome_to_taxid)
        assert_file_exists(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        with self.assertRaises(subprocess.CalledProcessError):
            _ = subprocess.run(cmd, shell=True, check=True)

    def test_wrong_outdir(self):
        yacht_output = PATH_TO_YACHT_OUTPUT
        genome_to_taxid = PATH_TO_GENOME_TO_TAXID
        outdir = os.path.join(os.path.dirname(__file__), 'testdata_nonexisting')

        assert_file_exists(yacht_output)
        assert_file_exists(genome_to_taxid)
        create_outdir(outdir)

        cmd = f"yacht convert --yacht_output {yacht_output} --sheet_name min_coverage0.2 --genome_to_taxid {genome_to_taxid} --outfile_prefix cami_result --outdir {outdir}"
        res = subprocess.run(cmd, shell=True, check=True)
        assert res.returncode == 0
        assert_file_exists(outdir)

        cleanup_outdir(outdir)

if __name__ == '__main__':
    unittest.main()