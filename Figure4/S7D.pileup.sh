coolpup.py data/DN3EV.cool data/DN3EV.uniq.bedpe --outname DN3EV.DN3EV.uniq.txt
coolpup.py data/DN3KO.cool data/DN3EV.uniq.bedpe --outname DN3KO.DN3EV.uniq.txt
coolpup.py data/DN3EV.cool data/DN3KO.uniq.bedpe --outname DN3EV.DN3KO.uniq.txt
coolpup.py data/DN3KO.cool data/DN3KO.uniq.bedpe --outname DN3KO.DN3KO.uniq.txt

plotpup.py DN3EV.*.txt DN3KO*.txt --row_names WT,KO --col_names WT,KO --output FigureS7D.pdf
