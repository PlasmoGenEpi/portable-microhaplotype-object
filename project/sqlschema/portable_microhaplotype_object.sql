-- # Class: "MarkerOfInterest" Description: "A specific genomic location of interest, e.g. drug resistance, or other phenotypical marker"
--     * Slot: id Description: 
--     * Slot: marker_location_id Description: the genomic location
-- # Class: "TargetInfo" Description: "Information about a specific target within a genome"
--     * Slot: id Description: 
--     * Slot: target_name Description: an identifier for this target
--     * Slot: gene_name Description: an identifier of the gene, if any, is being covered with this targeted
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
--     * Slot: insert_location_id Description: the intended genomic location of the insert of the amplicon (the location between the end of the forward primer and the beginning of the reverse primer)
-- # Class: "ReactionInfo" Description: "information on a panel of targeted amplicon primer pairs"
--     * Slot: id Description: 
--     * Slot: reaction_name Description: a name for this reaction
-- # Class: "PanelInfo" Description: "information on a panel of targeted amplicon primer pairs"
--     * Slot: id Description: 
--     * Slot: panel_name Description: a name for the panel
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "RepresentativeMicrohaplotype" Description: "the representative sequence for a microhaplotype, similar to a fast(a/q) format"
--     * Slot: id Description: 
--     * Slot: seq Description: the DNA sequence
--     * Slot: microhaplotype_name Description: an optional name for this microhaplotype
--     * Slot: quality Description: the ansi fastq per base quality score for this sequence, this is optional
--     * Slot: pseudocigar Description: the pseudocigar of the haplotype
-- # Class: "MaskingInfo" Description: "information about a subsegment of the sequence that should be masked"
--     * Slot: id Description: 
--     * Slot: seq_start Description: the start of the masking
--     * Slot: seq_segment_size Description: the size of the masking
--     * Slot: replacement_size Description: the size of replacement mask
-- # Class: "RepresentativeMicrohaplotypes" Description: "a collection of representative sequences for microhaplotypes for all targets"
--     * Slot: id Description: 
-- # Class: "RepresentativeMicrohaplotypesForTarget" Description: "a list of the representative sequence for a microhaplotypes, similar to a fast(a/q) format"
--     * Slot: target_id Description: the index into the target_info list
--     * Slot: RepresentativeMicrohaplotypes_id Description: Autocreated FK slot
-- # Class: "MicrohaplotypesDetected" Description: "the microhaplotypes detected in a targeted amplicon analysis"
--     * Slot: id Description: 
--     * Slot: bioinformatics_run_id Description: the index into bioinformatics_run_info list
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "GenomeInfo" Description: "information on a genome"
--     * Slot: id Description: 
--     * Slot: name Description: name of the genome
--     * Slot: genome_version Description: the genome version
--     * Slot: taxon_id Description: the NCBI taxonomy number
--     * Slot: url Description: a link to the where this genome file could be downloaded
--     * Slot: gff_url Description: a link to the where this genome's annotation file could be downloaded
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "GenomicLocation" Description: "information on the genomic location of specific sequence"
--     * Slot: id Description: 
--     * Slot: genome_id Description: the index to the genome in the targeted_genomes list that this location refers to 
--     * Slot: chrom Description: the chromosome name
--     * Slot: start Description: the start of the location, 0-based positioning
--     * Slot: end Description: the end of the location, 0-based positioning
--     * Slot: strand Description: which strand the location is, either + for plus strand or - for negative strand
--     * Slot: ref_seq Description: the reference sequence of this genomic location
-- # Class: "PrimerInfo" Description: "information on a primer sequence"
--     * Slot: id Description: 
--     * Slot: seq Description: the DNA sequence
--     * Slot: TargetInfo_id Description: Autocreated FK slot
--     * Slot: location_id Description: what the intended genomic location of the primer is
-- # Class: "MicrohaplotypesForSample" Description: "Microhaplotypes detected for a sample for all targets"
--     * Slot: id Description: 
--     * Slot: experiment_sample_id Description: the index into the experiment_info list
--     * Slot: MicrohaplotypesDetected_id Description: Autocreated FK slot
-- # Class: "MicrohaplotypeForTarget" Description: "Microhaplotype detected for a specific target"
--     * Slot: id Description: 
--     * Slot: mhap_id Description: the index for a microhaplotype for a target in the microhaplotypes_info list, e.g. microhaplotypes_info[mhaps_target_id][mhap_id]
--     * Slot: reads Description: the read count associated with this microhaplotype
--     * Slot: umis Description: the unique molecular identifier (umi) count associated with this microhaplotype
--     * Slot: MicrohaplotypesForTarget_id Description: Autocreated FK slot
-- # Class: "MicrohaplotypesForTarget" Description: "Microhaplotypes detected for a specific target"
--     * Slot: id Description: 
--     * Slot: mhaps_target_id Description: the index for a target in the microhaplotypes_info list
--     * Slot: MicrohaplotypesForSample_id Description: Autocreated FK slot
-- # Class: "BioinformaticsMethodInfo" Description: "the targeted amplicon bioinformatics pipeline"
--     * Slot: id Description: 
--     * Slot: bioinformatics_method_name Description: name of the collection of methods is called, e.g. pipeline 
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
--     * Slot: demultiplexing_method_id Description: the demultiplexing method used to separate raw reads from barcodes and primer targets
--     * Slot: denoising_method_id Description: the method used to de-noise and/or cluster the raw reads
-- # Class: "BioMethod" Description: "methodology description of a portion of a bioinformatics pipeline"
--     * Slot: id Description: 
--     * Slot: program_version Description: the version of generation method, should be in the format of v[MAJOR].[MINOR].[PATCH]
--     * Slot: program Description: name of the program used for this portion of the pipeline
--     * Slot: program_description Description: a short description of what this method does
--     * Slot: BioinformaticsMethodInfo_id Description: Autocreated FK slot
-- # Class: "ExperimentInfo" Description: "Information about a specific amplification and sequencing of a specimen"
--     * Slot: sequencing_info_id Description: the index into the sequencing_info list
--     * Slot: plate_name Description: a name of plate the specimen was in
--     * Slot: plate_row Description: the row the specimen was in
--     * Slot: plate_col Description: the column the specimen was in
--     * Slot: specimen_id Description: the index into the specimen_info list
--     * Slot: panel_id Description: the index into the panel_info list
--     * Slot: accession Description: ERA/SRA accession number for the sample if it was submitted
--     * Slot: experiment_sample_name Description: a unique identifier for this sequence/amplification run on a specimen_name
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "SequencingInfo" Description: "Information on sequencing info"
--     * Slot: id Description: 
--     * Slot: sequencing_info_name Description: a name of for the sequencing done, e.g. batch1
--     * Slot: seq_instrument Description: the sequencing instrument used to sequence the run, e.g. ILLUMINA, Illumina MiSeq
--     * Slot: seq_date Description: the date of sequencing, should be YYYY-MM or YYYY-MM-DD
--     * Slot: nucl_acid_ext Description: Link to a reference or kit that describes the recovery of nucleic acids from the sample
--     * Slot: nucl_acid_amp Description: Link to a reference or kit that describes the enzymatic amplification of nucleic acids
--     * Slot: nucl_acid_ext_date Description: the date of the nucleoacid extraction
--     * Slot: nucl_acid_amp_date Description: the date of the nucleoacid amplification
--     * Slot: pcr_cond Description: the method/conditions for PCR, List PCR cycles used to amplify the target
--     * Slot: lib_screen Description: Describe enrichment, screening, or normalization methods applied during amplification or library preparation, e.g. size selection 390bp, diluted to 1 ng DNA/sample
--     * Slot: lib_layout Description: Specify the configuration of reads, e.g. paired-end
--     * Slot: lib_kit Description: Name, version, and applicable cell or cycle numbers for the kit used to prepare libraries and load cells or chips for sequencing. If possible, include a part number, e.g. MiSeq Reagent Kit v3 (150-cycle), MS-102-3001
--     * Slot: seq_center Description: Name of facility where sequencing was performed (lab, core facility, or company)
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "ParasiteDensity" Description: "method and value of determined parasite density"
--     * Slot: id Description: 
--     * Slot: method Description: the method of how this density was obtained
--     * Slot: density Description: the density in microliters
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
-- # Class: "SpecimenInfo" Description: "Information on specimen info"
--     * Slot: plate_name Description: a name of plate the specimen was in
--     * Slot: plate_row Description: the row the specimen was in
--     * Slot: plate_col Description: the column the specimen was in
--     * Slot: specimen_name Description: an identifier for the specimen
--     * Slot: samp_taxon_id Description: the NCBI taxonomy number of the organism of interest
--     * Slot: individual_id Description: an identifier for the individual a specimen was collected from
--     * Slot: host_taxon_id Description: optional the NCBI taxonomy number of the host of the organism
--     * Slot: host_sex Description: if specimen is from a person, the sex listed for that person
--     * Slot: collection_date Description: the date of the sample collection, can be YYYY, YYYY-MM, or YYYY-MM-DD
--     * Slot: host_date_of_birth Description: if specimen is from a person, the date of birth of that person, can be YYYY, YYYY-MM, or YYYY-MM-DD
--     * Slot: collection_country Description: the name of country collected in, would be the same as admin level 0
--     * Slot: geo_admin1 Description: geographical admin level 1, the secondary large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin2 Description: geographical admin level 2, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin3 Description: geographical admin level 3, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: lat_lon Description: the latitude and longitude of the collection site of the specimen
--     * Slot: collector Description: the name of the primary person managing the specimen
--     * Slot: samp_store_loc Description: the sample store site, address or facility name
--     * Slot: samp_collect_device Description: the way the sample was collected, e.g. whole blood, dried blood spot, etc
--     * Slot: project_name Description: a name of the project under which the sample is organized
--     * Slot: sample_comments Description: any additional comments about the sample
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "BioinformaticsRunInfo" Description: "Information about the pipeline run that generated some of the microhaplotype detected and reads_by_stage"
--     * Slot: id Description: 
--     * Slot: bioinformatics_methods_id Description: the index into the bioinformatics_methods_info list
--     * Slot: run_date Description: the date when the run was done, should be YYYY-MM-DD
--     * Slot: bioinformatics_run_name Description: a name to for this run
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "PmoGenerationMethod" Description: "Information about how a PMO was generated"
--     * Slot: id Description: 
--     * Slot: program_version Description: the version of generation method, should be in the format of v[MAJOR].[MINOR].[PATCH]
--     * Slot: program_name Description: the name of the program
-- # Class: "PmoHeader" Description: "Information on the PMO file"
--     * Slot: id Description: 
--     * Slot: pmo_version Description: the version of the PMO file, should be in the format of v[MAJOR].[MINOR].[PATCH]
--     * Slot: creation_date Description: the date of when the PMO file was created or modified, should be YYYY-MM-DD
--     * Slot: generation_method_id Description: the generation method to create this PMO 
-- # Class: "StageReadCounts" Description: "Information on the reads counts at several stages"
--     * Slot: id Description: 
--     * Slot: read_count Description: the read counts
--     * Slot: stage Description: the stage of the pipeline, e.g. demultiplexed, denoised, etc
--     * Slot: ReadCountsByStageForTarget_id Description: Autocreated FK slot
-- # Class: "ReadCountsByStageForTarget" Description: "Information on the reads counts at several stages of a pipeline for a target"
--     * Slot: id Description: 
--     * Slot: target_id Description: the index into the target_info list
--     * Slot: ReadCountsByStageForExperimentalSample_id Description: Autocreated FK slot
-- # Class: "ReadCountsByStageForExperimentalSample" Description: "Information on the reads counts at several stages of a pipeline for a experimental_sample"
--     * Slot: id Description: 
--     * Slot: experiment_sample_id Description: the index into the experiment_info list
--     * Slot: total_raw_count Description: the raw counts off the sequencing machine that a sample began with
--     * Slot: ReadCountsByStage_id Description: Autocreated FK slot
-- # Class: "ReadCountsByStage" Description: "Information on the reads counts at several stages of a pipeline"
--     * Slot: id Description: 
--     * Slot: bioinformatics_run_id Description: the index into bioinformatics_run_info list
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "PortableMicrohaplotypeObject" Description: "Information on final results from a targeted amplicon analysis"
--     * Slot: id Description: 
--     * Slot: microhaplotypes_info_id Description: a list of the information on the representative microhaplotypes
--     * Slot: pmo_header_id Description: the PMO information for this file including version etc
-- # Class: "MarkerOfInterest_associations" Description: ""
--     * Slot: MarkerOfInterest_id Description: Autocreated FK slot
--     * Slot: associations Description: a list of associations with this marker, e.g. SP resistance, etc
-- # Class: "TargetInfo_target_attributes" Description: ""
--     * Slot: TargetInfo_id Description: Autocreated FK slot
--     * Slot: target_attributes Description: a list of classification type for the primer target
-- # Class: "ReactionInfo_panel_targets" Description: ""
--     * Slot: ReactionInfo_id Description: Autocreated FK slot
--     * Slot: panel_targets Description: a list of the target indexes in the target_info list
-- # Class: "PanelInfo_reactions" Description: ""
--     * Slot: PanelInfo_id Description: Autocreated FK slot
--     * Slot: reactions_id Description: a list of 1 or more reactions that this panel contains, each reactions list the targets that were amplified in that reaction, e.g. pool1, pool2
-- # Class: "RepresentativeMicrohaplotype_masking" Description: ""
--     * Slot: RepresentativeMicrohaplotype_id Description: Autocreated FK slot
--     * Slot: masking_id Description: masking info for the sequence
-- # Class: "RepresentativeMicrohaplotype_alt_annotations" Description: ""
--     * Slot: RepresentativeMicrohaplotype_id Description: Autocreated FK slot
--     * Slot: alt_annotations Description: a list of additional annotations associated with this microhaplotype, e.g. wildtype, amino acid changes etc
-- # Class: "RepresentativeMicrohaplotypesForTarget_microhaplotypes" Description: ""
--     * Slot: RepresentativeMicrohaplotypesForTarget_target_id Description: Autocreated FK slot
--     * Slot: microhaplotypes_id Description: a list of the microhaplotypes detected for a target 
-- # Class: "GenomeInfo_chromosomes" Description: ""
--     * Slot: GenomeInfo_id Description: Autocreated FK slot
--     * Slot: chromosomes Description: a list of chromosomes found within this genome
-- # Class: "BioMethod_additional_argument" Description: ""
--     * Slot: BioMethod_id Description: Autocreated FK slot
--     * Slot: additional_argument Description: any additional arguments that differ from the default
-- # Class: "SpecimenInfo_alternate_identifiers" Description: ""
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
--     * Slot: alternate_identifiers Description: a list of optional alternative names for the samples

CREATE TABLE "ReactionInfo" (
	id INTEGER NOT NULL, 
	reaction_name TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "RepresentativeMicrohaplotype" (
	id INTEGER NOT NULL, 
	seq TEXT NOT NULL, 
	microhaplotype_name TEXT, 
	quality TEXT, 
	pseudocigar TEXT, 
	PRIMARY KEY (id)
);
CREATE TABLE "MaskingInfo" (
	id INTEGER NOT NULL, 
	seq_start INTEGER NOT NULL, 
	seq_segment_size INTEGER NOT NULL, 
	replacement_size INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "RepresentativeMicrohaplotypes" (
	id INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "GenomicLocation" (
	id INTEGER NOT NULL, 
	genome_id INTEGER NOT NULL, 
	chrom TEXT NOT NULL, 
	start INTEGER NOT NULL, 
	"end" INTEGER NOT NULL, 
	strand TEXT, 
	ref_seq TEXT, 
	PRIMARY KEY (id)
);
CREATE TABLE "BioinformaticsMethodInfo" (
	id INTEGER NOT NULL, 
	bioinformatics_method_name TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	demultiplexing_method_id INTEGER NOT NULL, 
	denoising_method_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id), 
	FOREIGN KEY(demultiplexing_method_id) REFERENCES "BioMethod" (id), 
	FOREIGN KEY(denoising_method_id) REFERENCES "BioMethod" (id)
);
CREATE TABLE "BioMethod" (
	id INTEGER NOT NULL, 
	program_version TEXT NOT NULL, 
	program TEXT NOT NULL, 
	program_description TEXT, 
	"BioinformaticsMethodInfo_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("BioinformaticsMethodInfo_id") REFERENCES "BioinformaticsMethodInfo" (id)
);
CREATE TABLE "PmoGenerationMethod" (
	id INTEGER NOT NULL, 
	program_version TEXT NOT NULL, 
	program_name TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "MarkerOfInterest" (
	id INTEGER NOT NULL, 
	marker_location_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(marker_location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "RepresentativeMicrohaplotypesForTarget" (
	target_id INTEGER NOT NULL, 
	"RepresentativeMicrohaplotypes_id" INTEGER, 
	PRIMARY KEY (target_id), 
	FOREIGN KEY("RepresentativeMicrohaplotypes_id") REFERENCES "RepresentativeMicrohaplotypes" (id)
);
CREATE TABLE "PmoHeader" (
	id INTEGER NOT NULL, 
	pmo_version TEXT NOT NULL, 
	creation_date TEXT, 
	generation_method_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(generation_method_id) REFERENCES "PmoGenerationMethod" (id)
);
CREATE TABLE "ReactionInfo_panel_targets" (
	"ReactionInfo_id" INTEGER, 
	panel_targets INTEGER NOT NULL, 
	PRIMARY KEY ("ReactionInfo_id", panel_targets), 
	FOREIGN KEY("ReactionInfo_id") REFERENCES "ReactionInfo" (id)
);
CREATE TABLE "RepresentativeMicrohaplotype_masking" (
	"RepresentativeMicrohaplotype_id" INTEGER, 
	masking_id INTEGER, 
	PRIMARY KEY ("RepresentativeMicrohaplotype_id", masking_id), 
	FOREIGN KEY("RepresentativeMicrohaplotype_id") REFERENCES "RepresentativeMicrohaplotype" (id), 
	FOREIGN KEY(masking_id) REFERENCES "MaskingInfo" (id)
);
CREATE TABLE "RepresentativeMicrohaplotype_alt_annotations" (
	"RepresentativeMicrohaplotype_id" INTEGER, 
	alt_annotations TEXT, 
	PRIMARY KEY ("RepresentativeMicrohaplotype_id", alt_annotations), 
	FOREIGN KEY("RepresentativeMicrohaplotype_id") REFERENCES "RepresentativeMicrohaplotype" (id)
);
CREATE TABLE "BioMethod_additional_argument" (
	"BioMethod_id" INTEGER, 
	additional_argument TEXT, 
	PRIMARY KEY ("BioMethod_id", additional_argument), 
	FOREIGN KEY("BioMethod_id") REFERENCES "BioMethod" (id)
);
CREATE TABLE "PortableMicrohaplotypeObject" (
	id INTEGER NOT NULL, 
	microhaplotypes_info_id INTEGER NOT NULL, 
	pmo_header_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(microhaplotypes_info_id) REFERENCES "RepresentativeMicrohaplotypes" (id), 
	FOREIGN KEY(pmo_header_id) REFERENCES "PmoHeader" (id)
);
CREATE TABLE "MarkerOfInterest_associations" (
	"MarkerOfInterest_id" INTEGER, 
	associations TEXT, 
	PRIMARY KEY ("MarkerOfInterest_id", associations), 
	FOREIGN KEY("MarkerOfInterest_id") REFERENCES "MarkerOfInterest" (id)
);
CREATE TABLE "RepresentativeMicrohaplotypesForTarget_microhaplotypes" (
	"RepresentativeMicrohaplotypesForTarget_target_id" INTEGER, 
	microhaplotypes_id INTEGER NOT NULL, 
	PRIMARY KEY ("RepresentativeMicrohaplotypesForTarget_target_id", microhaplotypes_id), 
	FOREIGN KEY("RepresentativeMicrohaplotypesForTarget_target_id") REFERENCES "RepresentativeMicrohaplotypesForTarget" (target_id), 
	FOREIGN KEY(microhaplotypes_id) REFERENCES "RepresentativeMicrohaplotype" (id)
);
CREATE TABLE "TargetInfo" (
	id INTEGER NOT NULL, 
	target_name TEXT NOT NULL, 
	gene_name TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	insert_location_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id), 
	FOREIGN KEY(insert_location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "PanelInfo" (
	id INTEGER NOT NULL, 
	panel_name TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "MicrohaplotypesDetected" (
	id INTEGER NOT NULL, 
	bioinformatics_run_id INTEGER NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "GenomeInfo" (
	id INTEGER NOT NULL, 
	name TEXT NOT NULL, 
	genome_version TEXT NOT NULL, 
	taxon_id INTEGER NOT NULL, 
	url TEXT NOT NULL, 
	gff_url TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "ExperimentInfo" (
	sequencing_info_id INTEGER NOT NULL, 
	plate_name TEXT, 
	plate_row TEXT, 
	plate_col INTEGER, 
	specimen_id INTEGER NOT NULL, 
	panel_id INTEGER NOT NULL, 
	accession TEXT, 
	experiment_sample_name TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (experiment_sample_name), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "SequencingInfo" (
	id INTEGER NOT NULL, 
	sequencing_info_name TEXT NOT NULL, 
	seq_instrument TEXT NOT NULL, 
	seq_date TEXT NOT NULL, 
	nucl_acid_ext TEXT, 
	nucl_acid_amp TEXT, 
	nucl_acid_ext_date TEXT, 
	nucl_acid_amp_date TEXT, 
	pcr_cond TEXT, 
	lib_screen TEXT, 
	lib_layout TEXT, 
	lib_kit TEXT, 
	seq_center TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "SpecimenInfo" (
	plate_name TEXT, 
	plate_row TEXT, 
	plate_col INTEGER, 
	specimen_name TEXT NOT NULL, 
	samp_taxon_id INTEGER NOT NULL, 
	individual_id INTEGER, 
	host_taxon_id INTEGER, 
	host_sex TEXT, 
	collection_date TEXT NOT NULL, 
	host_date_of_birth TEXT, 
	collection_country TEXT NOT NULL, 
	geo_admin1 TEXT, 
	geo_admin2 TEXT, 
	geo_admin3 TEXT, 
	lat_lon TEXT, 
	collector TEXT NOT NULL, 
	samp_store_loc TEXT NOT NULL, 
	samp_collect_device TEXT NOT NULL, 
	project_name TEXT NOT NULL, 
	sample_comments TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (specimen_name), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "BioinformaticsRunInfo" (
	id INTEGER NOT NULL, 
	bioinformatics_methods_id INTEGER NOT NULL, 
	run_date TEXT NOT NULL, 
	bioinformatics_run_name TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "ReadCountsByStage" (
	id INTEGER NOT NULL, 
	bioinformatics_run_id INTEGER NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "PrimerInfo" (
	id INTEGER NOT NULL, 
	seq TEXT NOT NULL, 
	"TargetInfo_id" INTEGER, 
	location_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("TargetInfo_id") REFERENCES "TargetInfo" (id), 
	FOREIGN KEY(location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "MicrohaplotypesForSample" (
	id INTEGER NOT NULL, 
	experiment_sample_id INTEGER NOT NULL, 
	"MicrohaplotypesDetected_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("MicrohaplotypesDetected_id") REFERENCES "MicrohaplotypesDetected" (id)
);
CREATE TABLE "ParasiteDensity" (
	id INTEGER NOT NULL, 
	method TEXT NOT NULL, 
	density FLOAT NOT NULL, 
	"SpecimenInfo_specimen_name" TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "ReadCountsByStageForExperimentalSample" (
	id INTEGER NOT NULL, 
	experiment_sample_id INTEGER NOT NULL, 
	total_raw_count INTEGER NOT NULL, 
	"ReadCountsByStage_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("ReadCountsByStage_id") REFERENCES "ReadCountsByStage" (id)
);
CREATE TABLE "TargetInfo_target_attributes" (
	"TargetInfo_id" INTEGER, 
	target_attributes TEXT, 
	PRIMARY KEY ("TargetInfo_id", target_attributes), 
	FOREIGN KEY("TargetInfo_id") REFERENCES "TargetInfo" (id)
);
CREATE TABLE "PanelInfo_reactions" (
	"PanelInfo_id" INTEGER, 
	reactions_id INTEGER NOT NULL, 
	PRIMARY KEY ("PanelInfo_id", reactions_id), 
	FOREIGN KEY("PanelInfo_id") REFERENCES "PanelInfo" (id), 
	FOREIGN KEY(reactions_id) REFERENCES "ReactionInfo" (id)
);
CREATE TABLE "GenomeInfo_chromosomes" (
	"GenomeInfo_id" INTEGER, 
	chromosomes TEXT, 
	PRIMARY KEY ("GenomeInfo_id", chromosomes), 
	FOREIGN KEY("GenomeInfo_id") REFERENCES "GenomeInfo" (id)
);
CREATE TABLE "SpecimenInfo_alternate_identifiers" (
	"SpecimenInfo_specimen_name" TEXT, 
	alternate_identifiers TEXT, 
	PRIMARY KEY ("SpecimenInfo_specimen_name", alternate_identifiers), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "MicrohaplotypesForTarget" (
	id INTEGER NOT NULL, 
	mhaps_target_id INTEGER NOT NULL, 
	"MicrohaplotypesForSample_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("MicrohaplotypesForSample_id") REFERENCES "MicrohaplotypesForSample" (id)
);
CREATE TABLE "ReadCountsByStageForTarget" (
	id INTEGER NOT NULL, 
	target_id INTEGER NOT NULL, 
	"ReadCountsByStageForExperimentalSample_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("ReadCountsByStageForExperimentalSample_id") REFERENCES "ReadCountsByStageForExperimentalSample" (id)
);
CREATE TABLE "MicrohaplotypeForTarget" (
	id INTEGER NOT NULL, 
	mhap_id INTEGER NOT NULL, 
	reads INTEGER NOT NULL, 
	umis INTEGER, 
	"MicrohaplotypesForTarget_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("MicrohaplotypesForTarget_id") REFERENCES "MicrohaplotypesForTarget" (id)
);
CREATE TABLE "StageReadCounts" (
	id INTEGER NOT NULL, 
	read_count INTEGER NOT NULL, 
	stage TEXT NOT NULL, 
	"ReadCountsByStageForTarget_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("ReadCountsByStageForTarget_id") REFERENCES "ReadCountsByStageForTarget" (id)
);