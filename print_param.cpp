#include "print_param.h"

print_param::print_param(string file_Name)
{
    cout << "Starting up Parameter Print" << endl
         << endl;
    this->file_Name = file_Name;
}

void print_param::ingress()
{
    if (filesystem::exists(file_Name) == 0)
    {
        cout << "Sample parameter file being printed: " << file_Name << endl;

        fstream output;
        output.open(file_Name, ios::out);

        output << "{\n";

        output << "\t# Cuda device\n"
               << "\t\"CUDA Device ID\":0,\n\n";

        output << "\t# IO directory listings and universal file paths\n"
               << "\t#RECOMMENDS THE CREATION OF SEPERATE PARAMTER FILES WITH SEPERATE OUTPUT AND INTERMEDIATE PATHS FOR INDIVIDUAL PROJECTS\n "
               << "\t\"Input path\":\"results\",\n"
               << "\t\"Output path\":\"mk\",\n"
               << "\t\"Intermediate path\":\"intermediate\",\n\n";

        output << "\t# VCF sample details\n"
               << "\t\"Ploidy\":2,\n\n";

        output << "\t# Protein information\n"
               << "\t# DNA INTERPRETATION SHOULD BE USED FOR START AND STOP CODONS\n"
               << "\t\"Start codon(s)\":\"ATG\",\n"
               << "\t\"Stop codon(s)\":\"TAA,TAG,TGA\",\n"
               << "\t\"Genetic code\":\"A|GCT,GCC,GCA,GCG;R|CGT,CGC,CGA,CGG,AGA,AGG;N|AAT,AAC;D|GAT,GAC;B|AAT,AAC,GAT,GAC;C|TGT,TGC;Q|CAA,CAG;E|GAA,GAG;Z|CAA,CAG,GAA,GAG;G|GGT,GGC,GGA,GGG;H|CAT,CAC;M|ATG;I|ATT,ATC,ATA;L|CTT,CTC,CTA,CTG,TTA,TTG;K|AAA,AAG;F|TTT,TTC;P|CCT,CCC,CCA,CCG;S|TCT,TCC,TCA,TCG,AGT,AGC;T|ACT,ACC,ACA,ACG;W|TGG;Y|TAT,TAC;V|GTT,GTC,GTA,GTG;X|TAA,TGA,TAG\",\n\n";

        output << "\t# Any tab-deliminated text based formats such as but not limited to (*.txt, *.csv)\n"
               << "\t\"Universal gene list\":\"genelist.txt\",\n\n";

        output << "\t# VCF split parameters\n"
               << "\t\"Population file path\":\"/mnt/d/Deshan/Books/population/sample_population_codes.tsv\",\n"
               << "\t\"Reference allele count\":1,\n"
               << "\t\"Alternate allele count\":1,\n"
               << "\t\"SNP count per file\":10000,\n\n";

        output << "\t# FASTA split parameters\n"
               << "\t# REQUIRES SEPERATE INPUT, SINCE IT IS A SINGULAR INPUT\n"
               << "\t# IF ALL IS SET AS SEQUENCE THEN THE WHOLE FILE IS SPLIT. ELSE STATE THE SEQUENCE ID OF THE SEQUENCE THAT NEEDS TO BE SEPERATED\n"
               << "\t\"Sequence\":\"All\",\n"
               << "\t\"Raw FASTA file\":\"/mnt/d/Deshan/Books/University of Calgary/Experiments/Neutrality_Linux/reference_merge/One.fasta\",\n\n";

        output << "\t# FASTA merge parameters\n"
               << "\t# REQUIRES SEPERATE INPUT AND OUTPUT LOCATION\n"
               << "\t# ENSURE THE FASTA FILES HAVE THE APPROPRIATE EXTENSION: .fasta, .fna, .ffn, .faa, .frn, .fa\n"
               << "\t\"FASTA files folder\":\"/mnt/d/Deshan/Books/University of Calgary/Experiments/Neutrality_Linux/reference_Split/3\",\n"
               << "\t\"Merge FASTA path\" :\"/mnt/d/Deshan/Books/University of Calgary/Experiments/Neutrality_Linux/reference_Split/Three.fasta\",\n\n";

        output << "\t# All Neutrality tests\n"
                  "\t# Calculates all three neutrality tests(Tajima 's D, Fay and Wu' s and Fu and Li's) at once\n"
                  "\t# File created would have the extension *.nt (tab deliminated text file). File Name: CountryName_GeneListFileName.nt"
                  "\t\"Neutrality gene list\" : \"universal\",\n\n";

        output << "\t# Tajima's D\n"
               << "\t# File created would have the extension *.td (tab deliminated text file). File Name: CountryName_GeneListFileName.td\n"
               << "\t\"Tajima gene list\":\"universal\",\n\n";

        output << "\t# Fu and Li\n"
               << "\t# File created would have the extension *.fl (tab deliminated text file). File Name: CountryName_GeneListFileName.fl\n"
               << "\t\"Fu and Li gene list\":\"universal\",\n\n";

        output << "\t# Fay and Wu\n"
               << "\t# File created would have the extension *.fw (tab deliminated text file). File Name: CountryName_GeneListFileName.fw\n"
               << "\t\"Fay and Wu gene list\":\"universal\",\n\n";

        output << "\t# McDonald–Kreitman\n"
               << "\t# REFERENCE GENOME FILE SHOULD BE THE SAME AS THAT USED TO GENERATE THE VCF FILE\n"
               << "\t# ALIGNMENT FILE SHOULD BE A PAIRWISE ALIGNMENT OF THE REFERENCE GENOME TO THE OUTGROUP GENOME. SHOULD BE IN .maf FORMAT\n"
               << "\t# File created would have the extension *.mk (tab deliminated text file). File Name: CountryName_GeneListFileName.mk\n"
               << "\t\"Reference genome\":\"reference/Human_1.fasta\",\n"
               << "\t\"Alignment file\":\"reference/results.maf\",\n"
               << "\t\"McDonald–Kreitman gene list\":\"universal\",\n\n";

        output << "\t# Fst or Fixation Index\n"
                  "\t# Indexed population file should be a THREE column tab delimiated *.txt file with headings: Sample_name	population_ID	Super_population\n"
                  "\t# Ensure population_ID values are unique to each super population and do not overlap\n"
                  "\t\"Population index file path\":\"fst_pop/population.txt\",\n"
                  "\t\"Fst gene list\":\"universal\",\n"
                  "\t\"Population ID\":\"GWD,ACB,PJL\",\n\n";

        output << "\t# EHH or Extended Haplotype Homozygosity\n"
               << "\t# Define the extended haplotype region specification mode by either using a seperate gene list file (FILE mode) or state the standard displacement value (FIXED mode) in bases\n"
               << "\t\"Range mode\":\"FILE\",\n"
               << "\t\"EHH FILE path\":\"genelist_3.txt\",\n"
               << "\t\"FIXED mode\":\"+1000\"\n";

        output
            << "}";

        output.close();

        cout << "Sample parameter has been printed: " << file_Name << endl;
    }
    else
    {
        cout << "ERROR: PARAMETER FILE " << file_Name << " ALREADY EXISTS." << endl;
    }
}