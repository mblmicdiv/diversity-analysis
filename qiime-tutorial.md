# Microbial Diversity QIIME Lab

For this lab, we will be using QIIME (pronounced “chime”) for “Quantitative Insights Into Microbial Ecology”, an open source software package for comparison and analysis of microbial communities. The data used for QIIME is primarily based on high-throughput sequencing data (such as ssu rRNA diversity) but QIIME can also be used to analyze shotgun metagenomic data. QIIME has been cited by over 2,500 peer-reviewed journals since its [publication](http://www.nature.com/nmeth/journal/v7/n5/full/nmeth.f.303.html publication) in 2010. We will be running QIIME scripts using [MacQIIME](http://www.wernerlab.org/software/macqiime). It can be run this way on a Mac or as a Virtual Box on a PC. [QIIME installation](http://qiime.org/install/index.html)

It also can be run on an Amazon instance, which you will need to do if you have many or large datasets. [Information on running on Amazon instances](http://qiime.org/tutorials/working_with_aws.html)

From the microbiome samples, we will first call OTUs, assign taxonomy to the OTUs, and construct phylogenetic trees from representative sequences of OTUs. We will also compare the alpha diversity (within samples) and the beta diversity (between samples) by generating various charts and graphs that quantify diversity of the microbiome samples.

### Learning Goals for Microbiome Datasets using QIIME

- Understand steps in classifying and quantifying microbial diversity in environmental samples
- Classify various microbiome sequences into OTUs
- Compare and interpret rarefaction curves (alpha diversity estimates)
- Compare and interpret beta diversity of different treatments (cultivars and sampling regions)
- Develop hypotheses to explain potential differences in diversity between the various samples

### Quantitative Microbial Diversity Definitions

- **Microbiome**: all of the microbial organisms (bacteria, archaea, fungi, protists, viruses) that occur in a particular environment (or biome). Includes all of the genomes of these microbes and the interactions between these microbes.
- **Metagenomics**: genomic sequencing and analysis of all microbes from a particular environment without the need for cultivation; also termed “community” or “ecological” genomics.
- **OTU or operational taxonomic unit**: defined at 97% 16S rRNA sequence identity.
- **Alpha diversity**: description of microbial diversity ’‘within’’ a sample. Levels of alpha diversity are estimated using rarefaction plots.
- **Rarefaction curves**: reflect diversity within a sample based on abundance of various taxa within a community (OTU abundance vs. number of times sampled).
- **Beta diversity**: comparisons of diversity between a collection of samples. Metrics assess diversity differences between microbial communities. Output is a square “distance matrix” – visualized with PCOA or clustering. Weighted and unweighted unifrac – phylogenetic measures used to compare diversity between samples.

### Instructions for generating quantitative diversity analyses using QIIME (Virtual Box) ’

#### Installation

First- install [http://qiime.org/install/virtual_box.html QIIME VIrtual Box] on your machine (this may take > 1 hour - you can also use a machine in the lab) QIIME is installed on the Macs and PCs in the break room.

To get QIIME for your own computer:

Plug in your USB to a lab computer.
Transfer the folder c:\QIIME to your USB by dragging and dropping (or copy and pasting).
The folder should be ~12GB so the transfer will take 15 minutes or so.
Once the transfer completes, you should be able to double click the QIIME.vbox file to run the Virtual Machine
The ’‘‘password’’’ for the QIIME Virtual Box is qiime
QIIME is very RAM heavy. Therefore, some of the steps you’ll complete today will take a couple of minutes to run. Please be patient.
’‘‘Also’’’ - make sure you are logged into the shared drive:

“Finder” - > Go ->Location-> smb://nas.srv05-mblad.mbl.edu/md   
user: mblad/mm97924 pass: mm87924

—-
#### Get data

Attach to the 'md' server

Open the Terminal
Type:

    cd /Volumes/md
    ls
    cd MicDiv_2016
    ls
    cd Data_Repository
    cd Bioinformatics
    cd Day1_QIIME
    ls

The data we're using for this tutorial is in this directory.


### Classify and Summarize Diversity (Not covered in class: done separately)

#### Classify OTUs by making OTU table
 Operational taxonomical units (OTUs) are used to describe the various microbes in a sample. OTUs are defined as a cluster of reads with 97% 16S rRNA sequence identity. We will use QIIME to classify OTUs into an OTU table.

    pick_de_novo_otus.py -a -i seqs.fna -o otuParallel
or

     pick_de_novo_otus.py –i Seqs.fna –o otus

*Note:* this command takes ~over an hour to run and uses all RAM for the Virtual Box basically freezing the system until it’s completed. Thus we have run this step for you

#### View statistics of OTU table (SAVE)
View the statistics of the OTU table (a binary file) and save the output.

     biom summarize-table -i input-file.txt > otu_class.txt

From the OTU summary, look at how many OTUs correspond to each sample (“counts/sample detail”). What conclusions can you draw about the number of OTUs in these samples?

View and save PDF of barcharts (highest level = division). Find these in the wf_taxa_summary/taxa_summary_plots folder.

<!--
===Make OTU Heatmap===

Next, make a heat map of the OTUs per sample. The OTU table is visualized as a heat map where each row corresponds to an OTU and each column corresponds to a sample. The higher the relative abundance of an OTU in a sample, the more intense the color at the corresponding position in the heat map. OTUs are clustered by [http://en.wikipedia.org/wiki/UPGMA UPGMA hierarchical clustering]. QIIME indicates the biological classification by prefixing the level for example “p“ indicates phylum and ”g” indicates genus. For a refresher on biological classification, view this [http://en.wikipedia.org/wiki/Bacterial_taxonomy helpful wiki page].

time make_otu_heatmap.py -i otuParallel/otu_table.biom -o otuParallel/OTU_Heatmap.pdf
Precomputed one here: https://www.dropbox.com/s/15jp0su8i3lwzsv/OTU.heatmap.pdf?dl=0

Although, the resolution of the y-axis makes it difficult to read each OTU, it is still a valuable preliminary visualization. What types of information can you gain from this heat map? Are there any trends present at this stage with respect to the various samples?

===View and save PDF of the OTU Heatmap===

firefox otuParallel/OTU_Heatmap/*html
===Summarize community diversity by taxonomic composition (bar charts)===

Use the following command to generate box and area plots to describe our samples.

summarize_taxa_through_plots.py -i otuParallel/otu_table.biom -o otuParallel/wf_taxa_summary -m mappingFile.txt
-m provides the path to the mapping file with sample meta data

Double click the HTML files in the wf_taxa_summary/taxa_summary_plots folder. These should open in your browser. Take a look at the different plots and tables that are generated. If it’s hard to view the whole file, each of the plots are saved as PDFs in the charts folder within taxa_summary_plots. Which groups are the predominant phyla in the different samples? (Remember phyla is designated by “_p”.) Are there any predominant groups?

-->

### Determine the taxonomic diversity within each sample

#### Compute and compare alpha diversity – generate rarefaction curves

Alpha diversity tells us about the richness of species diversity within our sample.There are more than two dozen different established metrics to calculate the alpha diversity, we are using rarefaction. Feel free to read more details about other metrics [http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html here].

**Compute the alpha diversity by generating rarefaction curves.**


-m passes the mapping information such as treatment or cultivar-p passes the methods used to calculate the alpha diversity-t passes the path to the phylogenetic tree file

**Open** the **rarefaction_plots.html** and test different parameters to plot the metric tested against a category such as the sample site.

Is there an alpha diversity metric that estimates the completeness of our sample diversity differently than the other metrics? If so, which one is it? Does this metric estimate a higher or lower sample diversity compared to the other metrics?


    alpha_diversity.py -i otuParallel/otu_table.biom -m shannon,PD_whole_tree,chao1,observed_species -o otuParallel/wf_arare -t otuParallel/rep_set.tre

#### View and save PDF of rarefaction analyses

 Find this in the ’‘‘wf_arare/alpha_rarefaction_plots/rarefaction_plots.html.’’’ Then select PD_whole_tree and Treatment (or other combinations) *Compute beta diversity and compare beta diveristy by generating PCOA plots and UPGMA trees for Unifraq

#### Determine the Diversity Between Samples

The definition of beta diversity varies between ecologists. Here we will define beta diversity as the differentiation amongst habitats. To quantify beta diversity, QIIME calculates the pairwise dissimilarity between samples resulting in a distance matrix. For more information about other definitions/uses of beta diversity, see the [wikipedia page](http://en.wikipedia.org/wiki/Beta_diversity)

#### Compute the beta diversity and generate PCoA plots

    beta_diversity_through_plots.py -i otuParallel/otu_table.biom -o betaDiv -t otuParallel/rep_set.tre -m mapping.txt -e 100

*Note:* You may see an error about negative Eigenvalues, but the negative values are ~100x smaller than the positive values so we can ignore the warning.

-e sets the sequencing depth per sample
We are setting the sequencing depth of 289 is the minimum number of sequences among any sample. This will ensure there is not bias by some samples having higher read depth than others. This script returns a distance matrix and principal coordinate analysis (PCoA) plots. The dissimilarity between samples is measured by the UniFrac method which calculates the phylogenetic distance between sets of taxa. Weighted UniFrac (opposed to unweighted) accounts for the relative abundance of each taxa within the communities.

PCoA is also known as multidimensional scaling and you should be familiar with it by now from our earlier labs. If you’d like more information about the differences between PCoA and PCA, check out this [helpful blog post](http://occamstypewriter.org/boboh/2012/01/17/pca_and_pcoa_explained/) or [this course website](http://ordination.okstate.edu/overview.htm#Principal_coordinates_analysis).

**Open** the weighted and unweighted PCoA plots by double clicking the index.html in their respective folders. How does adjusting the PCoA plots for taxa abundance (weighted) affect the clustering and principal coordinates? On the colors tab, explore coloring by treatment, etc.

*Note:* These plots are use a lot of RAM. You could email yourself the zipped bf_bdiv_even289 folder and view them on the Windows machine, if the plots aren’t working well on the Virtual Box.

#### What are the significant correlations of particular samples?

The distance matrix generated for beta diversity can also be used to make UPGMA trees. UPGMA is a simple hierarchical clustering method and can be used to classify sampling units on the basis of pairwise similarities in dendrograms.

Now use the beta diversity to generate UPGMA trees.

    ￼upgma_cluster.py -i betaDiv/unweighted_unifrac_dm.txt -o unweighted.upgma.tre

    upgma_cluster.py –i betaDiv//weighted_unifrac_dm.txt –o weighted_upgma.tre

**Use phylodendron to visualize Unifrac trees** to view tree – open Newick formatted tree file in Firefox, and copy paste into [http://iubio.bio.indiana.edu/treeapp/treeprint-form.html] (http://iubio.bio.indiana.edu/treeapp/treeprint-form.html)

It’s most useful to view as a “phenogram” and save to PDF from within Phylodendron.

#### Analyses and PDFs you should have saved

OTU Heatmap

OTU taxa barcharts

Rarefaction analyses for alpha diversity

PCOA plot(s)

UPGMA Unifraq tree

#### Questions for the QIIME discussion

From the OTU table, how many OTUs correspond to each sample? What is the minimum/maximum of OTUs in all of the samples? Are there any significant differences between the samples? If so – why?

From the bar charts, which groups are the predominant phyla in the different samples? Are there any predominant groups unique to particular samples? Do you have any explanations for any of the observed differences between the different samples?

What types of information can you get from the “heat map”? How do you interpret the heat map with respect to the various microbiome samples?

From the rarefaction curves, have we fully sampled the diversity from the various sites? Why or why not?

From the UPGMA table – which sites cluster together (are more similar)? Which sites are different? How would you explain this pattern of diversity?

From the PCOA plots – are there any significant correlations of particular samples? If so, how would you explain these correlations?

Discuss some of the potential physiologies of the predominant groups of microbes correlated with the various microbome samples.

What do the differences in diversity tell you about the ecology of the different microbiome samples?
