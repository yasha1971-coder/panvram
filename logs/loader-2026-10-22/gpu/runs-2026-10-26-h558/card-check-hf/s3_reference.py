t0 = time.time()
REF = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz")
sh(f'curl -fsSL "{REF}" | gunzip > /content/cohort/reference.fa')
h = sh("sha256sum /content/cohort/reference.fa").split()[0]
assert h == "51ab7d3c2a375ecff53ee2a1be03e58d57cc1601c948a860c03a78d724f26e4c", h; print("reference sha256", h, flush=True)
stage("reference_download_verify", t0)
