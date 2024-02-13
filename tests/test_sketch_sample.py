import unittest
from yacht.sketch_sample import sketch_single_end, sketch_paired_end

def test_sketch_single_end(tmp_path):
    # Create a temporary file for testing
    infile = tmp_path / "test.fasta"
    infile.write_text(">seq1\nACGT")

    # Set up the arguments
    kmer = 31
    scaled = 1000
    outfile = tmp_path / "sketch.sig"

    # Call the function
    sketch_single_end(infile, kmer, scaled, outfile)

    # Assert that the output file exists
    assert outfile.exists()

    # Assert that the output file is not empty
    assert outfile.stat().st_size > 0


def test_sketch_paired_end(tmp_path):
    # Create temporary files for testing
    infile1 = tmp_path / "test1.fasta"
    infile1.write_text(">seq1\nACGT")
    infile2 = tmp_path / "test2.fasta"
    infile2.write_text(">seq2\nTGCA")

    # Set up the arguments
    kmer = 31
    scaled = 1000
    outfile = tmp_path / "sketch.sig"

    # Call the function
    sketch_paired_end(infile1, infile2, kmer, scaled, outfile)

    # Assert that the output file exists
    assert outfile.exists()

    # Assert that the output file is not empty
    assert outfile.stat().st_size > 0
