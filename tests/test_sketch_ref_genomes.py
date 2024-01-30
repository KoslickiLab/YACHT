import unittest
from yacht.sketch_ref_genomes import sketch_single_file, sketch_multiple_files

def test_sketch_single_file(tmp_path):
    # Create a temporary file for testing
    infile = tmp_path / "test.fasta"
    infile.write_text(">seq1\nACGT")

    # Set up the arguments
    kmer = 31
    scaled = 1000
    outfile = tmp_path / "sketch.sig"

    # Call the function
    sketch_single_file(infile, kmer, scaled, outfile)

    # Assert that the output file exists
    assert outfile.exists()

def test_sketch_multiple_files(tmp_path):
    # Create a temporary folder for testing
    folder_path = tmp_path / "test_folder"
    folder_path.mkdir()

    # Create some test files in the folder
    file1 = folder_path / "file1.fasta"
    file1.write_text(">seq1\nACGT")
    file2 = folder_path / "file2.fasta"
    file2.write_text(">seq2\nTGCA")

    # Set up the arguments
    kmer = 31
    scaled = 1000
    outfile = tmp_path / "sketch.sig"

    # Call the function
    sketch_multiple_files(folder_path, kmer, scaled, outfile)

    # Assert that the output file exists
    assert outfile.exists()
