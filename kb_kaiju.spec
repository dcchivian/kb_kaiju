/*
A KBase module: kb_kaiju
This SDK module is developed to wrap the open source package Kaiju
determines taxonomic structure for microbial communities from shotgun metagenomic sequence data

Module also utilizes Krona for visualization of results

References:
Kaiju Homepage: http://kaiju.binf.ku.dk/
Krona Homepage: https://github.com/marbl/Krona/wiki
Kaiju DBs from: http://kaiju.binf.ku.dk/server
Github repo for Kaiju: https://github.com/bioinformatics-centre/kaiju
Github repo for Krona: https://github.com/marbl/Krona

Kaiju paper: Menzel, P. et al. (2016) Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nat. Commun. 7:11257.
Krona paper: Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.

Kaiju License:
Copyright (c) 2015, 2016, 2017, 2018, 2019, 2020 Peter Menzel and Anders Krogh

Kaiju is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Kaiju is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the file LICENSE for more details.

You should have received a copy of the GNU General Public License along with the source code. If not, see http://www.gnu.org/licenses/.
*/

module kb_kaiju {
    /*
    A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int bool;

    /*
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string data_obj_name;
    typedef string data_obj_ref;


    /*
        kaiju command line:
	===================
	kaiju -t nodes.dmp -f kaiju_db_refseq.fmi -i inputfile.fastq  (for single-end libraries)
	kaiju -t nodes.dmp -f kaiju_db_refseq.fmi -i fwd_file.fastq -j rev_file.fastq  (for paired-end libraries)

	-o <output_file>
	-z <thread_count>                (e.g. 4)

	default run mode is MEM          (only consider exact matches)
	-a greedy -e <num_substitutions> (Greedy run mode permits <e> mismatches)
	-m <minimum_match_length>        (default: 11)
	-s <minimum_match_score>         (default: 65)

	Greedy mode:
	-E <e_value_cutoff>              (e.g. 0.05)

	-p                               (disable translation of input sequences to protein)
	-x                               (filter low-complexity regions with SEG - recommended)


	kaiju2table command line:
	=========================
	kaiju2table -t nodes.dmp -n names.dmp -r genus -o kaiju.out.summary kaiju.out

	-m <percent>                     (filter out members that are < percent of total reads)
	-u                               (filter out members that are < percent of all classified reads)
	-p                               (full taxon path instead of just taxon name)
  -l                               (taxon path containing only ranks specified by a comma-separated list)
  -c <integer > 0>                 (minimum required number of reads for the taxon (except viruses) to be reported)
  -e                               (expand viruses, which are always shown as full taxon path and read counts are not summarized in higher taxonomic levels)

	kaiju-addTaxonNames command line:
	===========================
	kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.names.out
	-u                               (omit unclassified reads)
	-p                               (print the full taxon path instead of just the taxon name)
	-r <class1>[,<class2>,...]       (print the path containing only to the specified ranks. e.g. "-r phylum,genus")

	kaiju2krona command lines:
	==========================
	kaiju2krona -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.out.krona
	ktImportText -o kaiju.out.html kaiju.out.krona

    */



    /* Kaiju App Input Params
    */
    typedef structure {
	workspace_name workspace_name;
	data_obj_ref   input_refs;
	list<string>   tax_levels;

	string         db_type;
	float          filter_percent;
	int            subsample_percent;
	int            subsample_replicates;
	int            subsample_seed;
	bool           seg_filter;
	int            min_match_length;
	bool           greedy_run_mode;
	int            greedy_allowed_mismatches;
	int            greedy_min_match_score;
	float          greedy_max_e_value;
	/*bool           filter_unclassified;*/  /* always filter unclassified */
	/*bool           full_tax_path;*/  /* will implement later */
	string         sort_taxa_by;
    } KaijuInputParams;


    /* Kaiju App Output
    */
    typedef structure {
        data_obj_name report_name;
        data_obj_ref  report_ref;
    } KaijuOutput;


    /* Kaiju Method
    */
    funcdef run_kaiju(KaijuInputParams params)
        returns (KaijuOutput) authentication required;
};
