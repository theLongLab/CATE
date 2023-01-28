#include "print_param.h"

print_param::print_param(string file_Name)
{
       /**
        * * Constructor Function
        * Assigns passed file location variable to the classes' private variable.
        **/

       cout << "Starting up Parameter Print" << endl
            << endl;
       this->file_Name = file_Name;
}

void print_param::ingress()
{
       /**
        * Execution function prints the paramater file to the user defined location in *.json format.
        **/
       // string original = this->file_Name;

       // for (size_t i = 1; i <= 22; i++)
       // {
       //        stringstream ss;
       //        ss << std::setw(2) << std::setfill('0') << i;
       //        string file_Name = original + ss.str() + ".json";

       /**
        * Ensures first and foremost that the entered file name does not exist in the said location.
        **/

       if (filesystem::exists(file_Name) == 0)
       {
              cout << "Sample parameter file being printed: " << file_Name << endl;

              fstream output;
              output.open(file_Name, ios::out);

              output << "{\n";

              output << "    # Cuda device\n"
                     << "    \"CUDA Device ID\":0,\n\n";

              output << "    # IO directory listings and universal file paths\n"
                     << "    # RECOMMENDS THE CREATION OF SEPARATE PARAMATER FILES WITH SEPARATE OUTPUT AND INTERMEDIATE PATHS FOR INDIVIDUAL PROJECTS\n"
                     << "    \"Input path\":\"VCF_indexed_folder\",\n"
                     << "    \"Output path\":\"Output_folder\",\n"
                     << "    \"Intermediate path\":\"Intermediate_folder\",\n\n";

              output << "    # Prometheus settings\n"
                     << "    # User is capable of determining the number of gene combinations and number of SNPs to be processed by CATE.\n"
                     << "    # When setting these parameters be mindful of system hardware availability.\n"
                     << "    # WARNING: Multi read CAN CAUSE A BOTTLENECK IF HARDWARE DOES NOT SUPPORT IT. NVMe DRIVES CONNECTED VIA PCIe BUS WILL SUPPORT MULTI READ. TYPICAL HDD DRIVES DO NOT.\n"
                     << "    \"Prometheus activate\":\"Yes\",\n"
                     << "    \"CPU cores\":25,\n"
                     << "    \"SNPs per time\":100000,\n"
                     << "    \"Number of genes\":1000,\n"
                     << "    \"Multi read\":\"Yes\",\n\n";

              output << "    # VCF sample details\n"
                     << "    \"Ploidy\":2,\n\n";

              output << "    # Protein information\n"
                     << "    # DNA INTERPRETATION SHOULD BE USED FOR START AND STOP CODONS\n"
                     << "    \"Start codon(s)\":\"ATG\",\n"
                     << "    \"Stop codon(s)\":\"TAA,TAG,TGA\",\n"
                     << "    \"Genetic code\":\"A|GCT,GCC,GCA,GCG;R|CGT,CGC,CGA,CGG,AGA,AGG;N|AAT,AAC;D|GAT,GAC;B|AAT,AAC,GAT,GAC;C|TGT,TGC;Q|CAA,CAG;E|GAA,GAG;Z|CAA,CAG,GAA,GAG;G|GGT,GGC,GGA,GGG;H|CAT,CAC;M|ATG;I|ATT,ATC,ATA;L|CTT,CTC,CTA,CTG,TTA,TTG;K|AAA,AAG;F|TTT,TTC;P|CCT,CCC,CCA,CCG;S|TCT,TCC,TCA,TCG,AGT,AGC;T|ACT,ACC,ACA,ACG;W|TGG;Y|TAT,TAC;V|GTT,GTC,GTA,GTG;X|TAA,TGA,TAG\",\n\n";

              output << "    # Neutrality and FST calculation mode\n"
                     << "    # Calculation mode can be FILE (to calculate the tests for predefined regions) or WINDOW\n"
                     << "    \"Calculation mode\":\"FILE\",\n\n";

              output << "    # WINDOW mode parameters\n"
                     << "    \"Window size\":10000,\n"
                     << "    \"Step size\":10000,\n\n";

              output << "    # FILE mode parameters\n"
                     << "    # Any tab-deliminated text based formats such as but not limited to (*.txt, *.csv)\n"
                     << "    \"Universal gene list\":\"sample.txt\n\n";

              output << "    # VCF split parameters\n"
                     << "    # The type of Split mode can be either CHR or CTSPLIT.\n"
                     << "    # The CHR mode will split the VCF file by chromosomes and extract the GT column only.\n"
                     << "    # The CTSPLIT mode will create CATE's proprietary file structure. For this the parent VCFs must only contain one chromosome's data and have only the GT column.\n"
                     << "    # Column numbers are non zero digits, i.e. Column numbers start with one.\n"
                     << "    \"Split mode\":\"CTSPLIT\"\n"
                     << "    \"CHR individual summary\":\"YES\"\n"
                     << "    \"Split cores\":10\n"
                     << "    \"Split SNPs per_time_CPU\":500000\n"
                     << "    \"Split SNPs per_time_GPU\":100000\n"
                     << "    \"Population file path\":\"sample_population_codes.tsv\",\n"
                     << "    \"Sample_ID Column number\":1,\n"
                     << "    \"Population_ID Column number\":6,\n"
                     << "    \"Reference allele count\":1,\n"
                     << "    \"Alternate allele count\":1,\n"
                     << "    \"SNP count per file\":10000,\n"
                     << "    \"MAF frequency\":\"0.05\"\n"
                     << "    \"Frequency logic\":\">\"\n\n";

              output << "    # FASTA split parameters\n"
                     << "    # REQUIRES SEPARATE INPUT, SINCE IT IS A SINGULAR INPUT\n"
                     << "    # IF ALL IS SET AS SEQUENCE THEN THE WHOLE FILE IS SPLIT. ELSE STATE THE SEQUENCE ID OF THE SEQUENCE THAT NEEDS TO BE SEPARATED\n"
                     << "    \"Sequence\":\"All\",\n"
                     << "    \"Raw FASTA file\":\"One.fasta\",\n\n";

              output << "    # FASTA merge parameters\n"
                     << "    # REQUIRES SEPARATE INPUT AND OUTPUT LOCATION\n"
                     << "    # ENSURE THE FASTA FILES HAVE THE APPROPRIATE EXTENSION: .fasta, .fna, .ffn, .faa, .frn, .fa\n"
                     << "    \"FASTA files folder\":\"reference_Split/\",\n"
                     << "    \"Merge FASTA path\":\"reference_Split/sample.fasta\",\n\n";

              output << "    # Extract genes parameters\n"
                     << "    # Extract gene FASTA sequences from REFERENCE FASTA file and outputs them as separate FASTA files\n"
                     << "    \"Reference genome ex\":\"reference/Human_1.fasta\",\n"
                     << "    \"Extract gene list\":\"universal\",\n\n";

              output << "    # GFF to Genes parameters\n"
                     << "    # REQUIRES SEPARATE INPUT, SINCE IT IS A SINGULAR INPUT\n"
                     << "    \"GFF file\":\"01.gff\",\n\n";

              output << "    # Haplotype INFO extract\n"
                     << "    # Extract haplotype information from VCFs\n"
                     << "    \"Reference genome hap\":\"reference/Human_1.fasta\",\n"
                     << "    \"Hap extract gene list\":\"universal\",\n\n";

              output << "    # All Neutrality tests\n"
                        "    # Calculates all three neutrality tests (Tajima's D, Fay and Wu's and Fu and Li's) at once\n"
                        "    # File created would have the extension *.nt (tab deliminated text file). File Name: CountryName_GeneListFileName.nt\n"
                        "    \"Neutrality gene list\":\"universal\",\n\n";

              output << "    # Tajima's D\n"
                     << "    # File created would have the extension *.td (tab deliminated text file). File Name: CountryName_GeneListFileName.td\n"
                     << "    \"Tajima gene list\":\"universal\",\n\n";

              output << "    # Fu and Li\n"
                     << "    # File created would have the extension *.fl (tab deliminated text file). File Name: CountryName_GeneListFileName.fl\n"
                     << "    \"Fu and Li gene list\":\"universal\",\n\n";

              output << "    # Fay and Wu\n"
                     << "    # File created would have the extension *.fw (tab deliminated text file). File Name: CountryName_GeneListFileName.fw\n"
                     << "    \"Fay and Wu gene list\":\"universal\",\n\n";

              output << "    # McDonald–Kreitman\n"
                     << "    # REFERENCE GENOME FILE SHOULD BE THE SAME AS THAT USED TO GENERATE THE VCF FILE\n"
                     << "    # ALIGNMENT FILE SHOULD BE A PAIRWISE ALIGNMENT OF THE REFERENCE GENOME TO THE OUTGROUP GENOME. SHOULD BE IN .maf FORMAT\n"
                     << "    # Alignment can be gene wide (GENE mode) or chromosome (CHROM mode) wide.\n"
                     << "    # GENE mode would require each gene to be separately aligned and the location of each alignment file to be added as a third column to the gene list file.\n"
                     << "    # File created would have the extension *.mk (tab deliminated text file). File Name: CountryName_GeneListFileName.mk\n"
                     << "    # CHROM mode would require the outgroup species chromosome's genome to be aligned with the entire chromosome of the query species reference genome.\n"
                     << "    # File created would have the extension *.mk (tab deliminated text file). File Name: CountryName_GeneListFileName.mk\n"
                     << "    \"Alignment mode\":\"GENE\",\n"
                     << "    \"ORF known\":\"Yes\",\n"
                     << "    \"Reference genome mk\":\"reference/Human_1.fasta\",\n"
                     << "    \"Alignment file\":\"reference/results.maf\",\n"
                     << "    \"McDonald–Kreitman gene list\":\"universal\",\n\n";

              output << "    # Fst or Fixation Index\n"
                        "    # Indexed population file should be a THREE column tab deliminated *.txt file with headings: Sample_name	population_ID	Super_population\n"
                        "    # Ensure population_ID values are unique to each super population and do not overlap\n"
                        "    \"Population index file path\":\"fst_pop/population.txt\",\n"
                        "    \"Fst gene list\":\"universal\",\n"
                        "    \"Population ID\":\"GWD,ACB,PJL\",\n\n";

              output << "    # EHH or Extended Haplotype Homozygosity\n"
                     << "    # Define the extended haplotype region specification mode by either using a separate gene list file (FILE mode) or state the standard displacement value (FIXED mode) in bases\n"
                     << "    # Degradation of EHH around a SNP can be calculated by SNP or BP mode.\n"
                     << "    \"Range mode\":\"FILE\",\n"
                     << "    \"EHH FILE path\":\"universal\",\n"
                     << "    \"FIXED mode\":\"+1000\"\n"
                     << "    \"SNP default count\":53\n"
                     << "    \"SNP BP displacement\":100000\n"
                     << "    \"EHH CPU cores\":10\n";

              output
                  << "}";

              output.close();

              cout << "Sample parameter has been printed: " << file_Name << endl;
       }
       else
       {
              /**
               *  ! Initialised if the parameter file already exists in the said location.
               **/

              cout << "ERROR: PARAMETER FILE " << file_Name << " ALREADY EXISTS." << endl;
       }
}