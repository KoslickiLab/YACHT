# Create the example sample data for a patient with respiratory symptoms seeks to find out the pathogen that is causing them these symptoms.

# Before moving on. Make sure reads needed to create sample dataset are available. Please reference create_reference_database.md

# Create samples that will be loaded to the 96-well tray

# Negative control, so just reads from a healthy lung
cat SRR2830253.fasta negative_control_well_11.fasta

# Positive control with H. influenzae
cat SRR25626360.fasta SRR2830253.fasta > positive_control_well_23.fasta

# Sample 1
cat SRR25626360.fasta SRR2830253.fasta SRR25626360.fasta > positive_control_well_64.fasta

# Sample 2
cat SRR24210460.fasta SRR2830253.fasta SRR25626360.fasta > sample_well_80.fasta

# I check one of my negative controls, which is a healthy lung example and we should not detect any bacteria here
# no contamination

# I check one of my positive controls for M. pneumonaie which should not have  H. influenzae
# no contamination

# I check one of my positive controls for H. influenzae which should not have M. pneumonaie
# contamination

# I check one of my samples for H. influenzae which should not have M. pneumonaie
# contamination

