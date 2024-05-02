import scanpy as sc

outpath = '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/Scrublet/Lai'

dirpaths = [
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR9897621_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR9897622_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR9897623_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR9897624_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR9897625_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR12539462_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR12539463_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR12603780_filtered.h5',
    '/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/CellBender/Lai/SRR14615558_filtered.h5'
]

for dirpath in dirpaths:
    id = dirpath.split('/')[10].split('_')[0]  # Define id here
    adata = sc.read_10x_h5(dirpath)  # Use dirpath directly here
    sc.external.pp.scrublet(adata, verbose=True)
    sc.write(f"{outpath}/{id}_scrub.h5", adata)