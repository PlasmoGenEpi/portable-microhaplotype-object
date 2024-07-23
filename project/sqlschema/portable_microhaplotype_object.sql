-- # Class: "TargetInfo" Description: "Information about a specific target within a genome"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
--     * Slot: gene_id Description: an identifier of the gene, if any, is being covered with this targeted
--     * Slot: PanelInfo_id Description: Autocreated FK slot
--     * Slot: insert_location_id Description: the intended genomic location of the insert of the amplicon (the location between the end of the forward primer and the beginning of the reverse primer)
-- # Class: "PanelInfo" Description: "information on a panel of targeted amplicon primer pairs"
--     * Slot: id Description: 
--     * Slot: panel_id Description: name of the panel
--     * Slot: target_genome_id Description: the info on the target reference genome for this panel
-- # Class: "RepresentativeMicrohaplotypeSequence" Description: "the representative sequence for a microhaplotype, similar to a fast(a/q) format"
--     * Slot: id Description: 
--     * Slot: seq Description: the DNA sequence
--     * Slot: microhaplotype_id Description: name of the microhaplotype, should be unique to this microhaplotype
--     * Slot: quality Description: the ansi fastq per base quality score for this sequence, this is optional
--     * Slot: pseudocigar Description: the pseudocigar of the haplotype
-- # Class: "MaskingInfo" Description: "information about a subsegment of the sequence that should be masked"
--     * Slot: id Description: 
--     * Slot: seq_start Description: the start of the masking
--     * Slot: seq_segment_size Description: the size of the masking
--     * Slot: replacement_size Description: the size of replacement mask
-- # Class: "RepresentativeMicrohaplotypeSequences" Description: "a collection of representative sequences for microhaplotypess for all targets"
--     * Slot: tar_amp_bioinformatics_info_id Description: a unique identifier for this targeted amplicon bioinformatics pipeline run
--     * Slot: PortableMicrohaplotypeObject_analysis_name Description: Autocreated FK slot
-- # Class: "RepresentativeMicrohaplotypesForTarget" Description: "a list of the representative sequence for a microhaplotypes, similar to a fast(a/q) format"
--     * Slot: target_id Description: name of the target
--     * Slot: RepresentativeMicrohaplotypeSequences_tar_amp_bioinformatics_info_id Description: Autocreated FK slot
-- # Class: "MicrohaplotypesDetected" Description: "the microhaplotypes detected in a targeted amplicon analysis"
--     * Slot: id Description: 
--     * Slot: tar_amp_bioinformatics_info_id Description: a unique identifier for this targeted amplicon bioinformatics pipeline run
-- # Class: "DemultiplexedExperimentSamples" Description: "a list of raw reads counts for each experiment sample for all targets within panel"
--     * Slot: id Description: 
--     * Slot: tar_amp_bioinformatics_info_id Description: a unique identifier for this targeted amplicon bioinformatics pipeline run
-- # Class: "DemultiplexedTargetsForExperimentSample" Description: "a list of raw reads for a experiment sample for all targets within panel"
--     * Slot: id Description: 
--     * Slot: experiment_sample_id Description: a unique identifier for this sequence/amplification run on a specimen
-- # Class: "DemultiplexedTargetForExperimentSample" Description: "the raw read count for a experiment sample for a target"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
--     * Slot: raw_read_count Description: the raw read counts extracted for a target for a experiment sample
-- # Class: "GenomeInfo" Description: "information on a genome"
--     * Slot: id Description: 
--     * Slot: name Description: name of the genome
--     * Slot: version Description: the genome version
--     * Slot: taxon_id Description: the NCBI taxonomy number
--     * Slot: url Description: a link to the where this genome file could be downloaded
--     * Slot: gff_url Description: a link to the where this genome's annotation file could be downloaded
-- # Class: "GenomicLocation" Description: "information on the genomic location of specific sequence"
--     * Slot: id Description: 
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
--     * Slot: experiment_sample_id Description: a unique identifier for this sequence/amplification run on a specimen
-- # Class: "MicrohaplotypeForTarget" Description: "Microhaplotype detected for a specific target"
--     * Slot: id Description: 
--     * Slot: microhaplotype_id Description: name of the microhaplotype, should be unique to this microhaplotype
--     * Slot: read_count Description: the read count associated with this microhaplotype
--     * Slot: umi_count Description: the unique molecular identifier (umi) count associated with this microhaplotype
--     * Slot: MicrohaplotypesForTarget_id Description: Autocreated FK slot
-- # Class: "MicrohaplotypesForTarget" Description: "Microhaplotypes detected for a specific target"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
-- # Class: "TarAmpBioinformaticsInfo" Description: "the targeted amplicon bioinformatics pipeline"
--     * Slot: tar_amp_bioinformatics_info_id Description: a unique identifier for this targeted amplicon bioinformatics pipeline run
--     * Slot: demultiplexing_method_id Description: the demultiplexing method used to separate raw reads from barcodes and primer targets
--     * Slot: denoising_method_id Description: the method used to de-noise and/or cluster the raw reads
--     * Slot: population_clustering_method_id Description: the method used to compare clustered/de-noised reads across samples for a target
-- # Class: "BioMethod" Description: "methodology description of a portion of a bioinformatics pipeline"
--     * Slot: id Description: 
--     * Slot: program Description: name of the program used for this portion of the pipeline
--     * Slot: purpose Description: the propose for this method
--     * Slot: program_version Description: versioning info for the program
--     * Slot: TarAmpBioinformaticsInfo_tar_amp_bioinformatics_info_id Description: Autocreated FK slot
-- # Class: "ExperimentInfo" Description: "Information about a specific amplification and sequencing of a specimen"
--     * Slot: sequencing_info_id Description: a unique identifier for this sequencing info
--     * Slot: plate_name Description: a name of plate the specimen was in
--     * Slot: plate_row Description: the row the specimen was in
--     * Slot: plate_col Description: the column the specimen was in
--     * Slot: specimen_id Description: the name of the specimen of a individual
--     * Slot: panel_id Description: name of the panel
--     * Slot: experiment_sample_id Description: a unique identifier for this sequence/amplification run on a specimen
--     * Slot: PortableMicrohaplotypeObject_analysis_name Description: Autocreated FK slot
-- # Class: "SequencingInfo" Description: "Information on sequencing info"
--     * Slot: sequencing_info_id Description: a unique identifier for this sequencing info
--     * Slot: seq_instrument Description: the sequencing instrument used to sequence the run, e.g. ILLUMINA, Illumina MiSeq
--     * Slot: seq_date Description: the date of sequencing, should be YYYY-MM or YYYY-MM-DD
--     * Slot: nucl_acid_ext Description: Link to a reference or kit that describes the recovery of nucleic acids from the sample
--     * Slot: nucl_acid_amp Description: Link to a reference or kit that describes the enzymatic amplification of nucleic acids,
--     * Slot: nucl_acid_ext_date Description: the date of the nucleoacide extraction
--     * Slot: nucl_acid_amp_date Description: the date of the nucleoacide amplification
--     * Slot: pcr_cond Description: the method/conditions for PCR, List PCR cycles used to amplify the target
--     * Slot: lib_screen Description: Describe enrichment, screening, or normalization methods applied during amplification or library preparation, e.g. size selection 390bp, diluted to 1 ng DNA/sample
--     * Slot: lib_layout Description: Specify the configuration of reads, e.g. paired-end
--     * Slot: lib_kit Description: Name, version, and applicable cell or cycle numbers for the kit used to prepare libraries and load cells or chips for sequencing. If possible, include a part number, e.g. MiSeq Reagent Kit v3 (150-cycle), MS-102-3001
--     * Slot: seq_center Description: Name of facility where sequencing was performed (lab, core facility, or company)
-- # Class: "SpecimenInfo" Description: "Information on specimen info"
--     * Slot: specimen_id Description: the name of the specimen of a individual
--     * Slot: plate_name Description: a name of plate the specimen was in
--     * Slot: plate_row Description: the row the specimen was in
--     * Slot: plate_col Description: the column the specimen was in
--     * Slot: samp_taxon_id Description: the NCBI taxonomy number of the organism of interest
--     * Slot: individual_id Description: an identifier for the individual a specimen was collected from
--     * Slot: host_taxon_id Description: optional the NCBI taxonomy number of the host of the organism
--     * Slot: parasite_density Description: the parasite density in microliters
--     * Slot: collection_date Description: the date of the sample collection
--     * Slot: collection_country Description: the name of country collected in, would be the same as admin level 0
--     * Slot: geo_admin1 Description: geographical admin level 1, the secondary large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin2 Description: geographical admin level 2, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin3 Description: geographical admin level 3, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: lat_lon Description: the latitude and longitude of the collection site of the specimen
--     * Slot: collector Description: the name of the primary person managing the specimen
--     * Slot: samp_store_loc Description: the sample store site, address or facility name
--     * Slot: samp_collect_device Description: the way the sample was collected, e.g. whole blood, dried blood spot, etc
--     * Slot: project_name Description: a name of the project under which the sample is organized
--     * Slot: accession Description: ERA/SRA accession number for the sample if it was submitted
--     * Slot: sample_comments Description: any additional comments about the sample
--     * Slot: PortableMicrohaplotypeObject_analysis_name Description: Autocreated FK slot
-- # Class: "PortableMicrohaplotypeObject" Description: "Information on final results from a targeted amplicon analysis"
--     * Slot: analysis_name Description: a name for the analysis detailed here in this experiment, can be a concatenation of names if combined more than one PMO
--     * Slot: sequencing_infos_sequencing_info_id Description: the sequencing info for this project
--     * Slot: panel_info_id Description: the info on the panel used for this project
--     * Slot: microhaplotypes_detected_id Description: the microhaplotypes detected in this projects
--     * Slot: target_demultiplexed_experiment_samples_id Description: the raw demultiplex target counts for each sample
--     * Slot: taramp_bioinformatics_infos_tar_amp_bioinformatics_info_id Description: the bioinformatics pipeline/methods used to generated the amplicon analysis for this project
--     * Slot: postprocessing_bioinformatics_infos_id Description: any additional methods that were applied to the processing of this file/analysis, this can be filtering, adding info etc
-- # Class: "RepresentativeMicrohaplotypeSequence_masking" Description: ""
--     * Slot: RepresentativeMicrohaplotypeSequence_id Description: Autocreated FK slot
--     * Slot: masking_id Description: masking info for the sequence
-- # Class: "RepresentativeMicrohaplotypeSequence_alt_annotations" Description: ""
--     * Slot: RepresentativeMicrohaplotypeSequence_id Description: Autocreated FK slot
--     * Slot: alt_annotations Description: a list of additional annotations associated with this microhaplotype, e.g. wildtype, amino acid changes etc
-- # Class: "RepresentativeMicrohaplotypesForTarget_seqs" Description: ""
--     * Slot: RepresentativeMicrohaplotypesForTarget_target_id Description: Autocreated FK slot
--     * Slot: seqs_id Description: a list of the microhaplotypes detected for a target 
-- # Class: "MicrohaplotypesDetected_experiment_samples" Description: ""
--     * Slot: MicrohaplotypesDetected_id Description: Autocreated FK slot
--     * Slot: experiment_samples_id Description: a list of the microhaplotypes detected for a sample for various targets 
-- # Class: "DemultiplexedExperimentSamples_demultiplexed_experiment_samples" Description: ""
--     * Slot: DemultiplexedExperimentSamples_id Description: Autocreated FK slot
--     * Slot: demultiplexed_experiment_samples_id Description: a list of the samples with the number of raw reads extracted 
-- # Class: "DemultiplexedTargetsForExperimentSample_demultiplexed_targets" Description: ""
--     * Slot: DemultiplexedTargetsForExperimentSample_id Description: Autocreated FK slot
--     * Slot: demultiplexed_targets_id Description: a list of the targets extracted for a sample 
-- # Class: "MicrohaplotypesForSample_target_results" Description: ""
--     * Slot: MicrohaplotypesForSample_id Description: Autocreated FK slot
--     * Slot: target_results_id Description: a list of the microhaplotypes detected for a list of targets
-- # Class: "BioMethod_additional_argument" Description: ""
--     * Slot: BioMethod_id Description: Autocreated FK slot
--     * Slot: additional_argument Description: any additional arguments that differ from the default
-- # Class: "SpecimenInfo_alternate_identifiers" Description: ""
--     * Slot: SpecimenInfo_specimen_id Description: Autocreated FK slot
--     * Slot: alternate_identifiers Description: a list of optional alternative names for the samples

CREATE TABLE "RepresentativeMicrohaplotypeSequence" (
	id INTEGER NOT NULL, 
	seq TEXT NOT NULL, 
	microhaplotype_id TEXT NOT NULL, 
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
CREATE TABLE "MicrohaplotypesDetected" (
	id INTEGER NOT NULL, 
	tar_amp_bioinformatics_info_id TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "DemultiplexedExperimentSamples" (
	id INTEGER NOT NULL, 
	tar_amp_bioinformatics_info_id TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "DemultiplexedTargetsForExperimentSample" (
	id INTEGER NOT NULL, 
	experiment_sample_id TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "DemultiplexedTargetForExperimentSample" (
	id INTEGER NOT NULL, 
	target_id TEXT NOT NULL, 
	raw_read_count FLOAT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "GenomeInfo" (
	id INTEGER NOT NULL, 
	name TEXT NOT NULL, 
	version TEXT NOT NULL, 
	taxon_id INTEGER NOT NULL, 
	url TEXT NOT NULL, 
	gff_url TEXT, 
	PRIMARY KEY (id)
);
CREATE TABLE "GenomicLocation" (
	id INTEGER NOT NULL, 
	chrom TEXT NOT NULL, 
	start INTEGER NOT NULL, 
	"end" INTEGER NOT NULL, 
	strand TEXT, 
	ref_seq TEXT, 
	PRIMARY KEY (id)
);
CREATE TABLE "MicrohaplotypesForSample" (
	id INTEGER NOT NULL, 
	experiment_sample_id TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "MicrohaplotypesForTarget" (
	id INTEGER NOT NULL, 
	target_id TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "TarAmpBioinformaticsInfo" (
	tar_amp_bioinformatics_info_id TEXT NOT NULL, 
	demultiplexing_method_id INTEGER NOT NULL, 
	denoising_method_id INTEGER NOT NULL, 
	population_clustering_method_id INTEGER NOT NULL, 
	PRIMARY KEY (tar_amp_bioinformatics_info_id), 
	FOREIGN KEY(demultiplexing_method_id) REFERENCES "BioMethod" (id), 
	FOREIGN KEY(denoising_method_id) REFERENCES "BioMethod" (id), 
	FOREIGN KEY(population_clustering_method_id) REFERENCES "BioMethod" (id)
);
CREATE TABLE "BioMethod" (
	id INTEGER NOT NULL, 
	program TEXT NOT NULL, 
	purpose TEXT NOT NULL, 
	program_version TEXT NOT NULL, 
	"TarAmpBioinformaticsInfo_tar_amp_bioinformatics_info_id" TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY("TarAmpBioinformaticsInfo_tar_amp_bioinformatics_info_id") REFERENCES "TarAmpBioinformaticsInfo" (tar_amp_bioinformatics_info_id)
);
CREATE TABLE "SequencingInfo" (
	sequencing_info_id TEXT NOT NULL, 
	seq_instrument TEXT NOT NULL, 
	seq_date TEXT NOT NULL, 
	nucl_acid_ext TEXT NOT NULL, 
	nucl_acid_amp TEXT NOT NULL, 
	nucl_acid_ext_date TEXT NOT NULL, 
	nucl_acid_amp_date TEXT NOT NULL, 
	pcr_cond TEXT NOT NULL, 
	lib_screen TEXT NOT NULL, 
	lib_layout TEXT NOT NULL, 
	lib_kit TEXT NOT NULL, 
	seq_center TEXT NOT NULL, 
	PRIMARY KEY (sequencing_info_id)
);
CREATE TABLE "PanelInfo" (
	id INTEGER NOT NULL, 
	panel_id TEXT NOT NULL, 
	target_genome_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(target_genome_id) REFERENCES "GenomeInfo" (id)
);
CREATE TABLE "MicrohaplotypeForTarget" (
	id INTEGER NOT NULL, 
	microhaplotype_id TEXT NOT NULL, 
	read_count FLOAT NOT NULL, 
	umi_count FLOAT, 
	"MicrohaplotypesForTarget_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("MicrohaplotypesForTarget_id") REFERENCES "MicrohaplotypesForTarget" (id)
);
CREATE TABLE "RepresentativeMicrohaplotypeSequence_masking" (
	"RepresentativeMicrohaplotypeSequence_id" INTEGER, 
	masking_id INTEGER, 
	PRIMARY KEY ("RepresentativeMicrohaplotypeSequence_id", masking_id), 
	FOREIGN KEY("RepresentativeMicrohaplotypeSequence_id") REFERENCES "RepresentativeMicrohaplotypeSequence" (id), 
	FOREIGN KEY(masking_id) REFERENCES "MaskingInfo" (id)
);
CREATE TABLE "RepresentativeMicrohaplotypeSequence_alt_annotations" (
	"RepresentativeMicrohaplotypeSequence_id" INTEGER, 
	alt_annotations TEXT, 
	PRIMARY KEY ("RepresentativeMicrohaplotypeSequence_id", alt_annotations), 
	FOREIGN KEY("RepresentativeMicrohaplotypeSequence_id") REFERENCES "RepresentativeMicrohaplotypeSequence" (id)
);
CREATE TABLE "MicrohaplotypesDetected_experiment_samples" (
	"MicrohaplotypesDetected_id" INTEGER, 
	experiment_samples_id INTEGER NOT NULL, 
	PRIMARY KEY ("MicrohaplotypesDetected_id", experiment_samples_id), 
	FOREIGN KEY("MicrohaplotypesDetected_id") REFERENCES "MicrohaplotypesDetected" (id), 
	FOREIGN KEY(experiment_samples_id) REFERENCES "MicrohaplotypesForSample" (id)
);
CREATE TABLE "DemultiplexedExperimentSamples_demultiplexed_experiment_samples" (
	"DemultiplexedExperimentSamples_id" INTEGER, 
	demultiplexed_experiment_samples_id INTEGER NOT NULL, 
	PRIMARY KEY ("DemultiplexedExperimentSamples_id", demultiplexed_experiment_samples_id), 
	FOREIGN KEY("DemultiplexedExperimentSamples_id") REFERENCES "DemultiplexedExperimentSamples" (id), 
	FOREIGN KEY(demultiplexed_experiment_samples_id) REFERENCES "DemultiplexedTargetsForExperimentSample" (id)
);
CREATE TABLE "DemultiplexedTargetsForExperimentSample_demultiplexed_targets" (
	"DemultiplexedTargetsForExperimentSample_id" INTEGER, 
	demultiplexed_targets_id INTEGER NOT NULL, 
	PRIMARY KEY ("DemultiplexedTargetsForExperimentSample_id", demultiplexed_targets_id), 
	FOREIGN KEY("DemultiplexedTargetsForExperimentSample_id") REFERENCES "DemultiplexedTargetsForExperimentSample" (id), 
	FOREIGN KEY(demultiplexed_targets_id) REFERENCES "DemultiplexedTargetForExperimentSample" (id)
);
CREATE TABLE "MicrohaplotypesForSample_target_results" (
	"MicrohaplotypesForSample_id" INTEGER, 
	target_results_id INTEGER NOT NULL, 
	PRIMARY KEY ("MicrohaplotypesForSample_id", target_results_id), 
	FOREIGN KEY("MicrohaplotypesForSample_id") REFERENCES "MicrohaplotypesForSample" (id), 
	FOREIGN KEY(target_results_id) REFERENCES "MicrohaplotypesForTarget" (id)
);
CREATE TABLE "BioMethod_additional_argument" (
	"BioMethod_id" INTEGER, 
	additional_argument TEXT, 
	PRIMARY KEY ("BioMethod_id", additional_argument), 
	FOREIGN KEY("BioMethod_id") REFERENCES "BioMethod" (id)
);
CREATE TABLE "TargetInfo" (
	id INTEGER NOT NULL, 
	target_id TEXT NOT NULL, 
	gene_id TEXT, 
	"PanelInfo_id" INTEGER, 
	insert_location_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PanelInfo_id") REFERENCES "PanelInfo" (id), 
	FOREIGN KEY(insert_location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "PortableMicrohaplotypeObject" (
	analysis_name TEXT NOT NULL, 
	sequencing_infos_sequencing_info_id TEXT NOT NULL, 
	panel_info_id INTEGER NOT NULL, 
	microhaplotypes_detected_id INTEGER NOT NULL, 
	target_demultiplexed_experiment_samples_id INTEGER, 
	taramp_bioinformatics_infos_tar_amp_bioinformatics_info_id TEXT NOT NULL, 
	postprocessing_bioinformatics_infos_id INTEGER, 
	PRIMARY KEY (analysis_name), 
	FOREIGN KEY(sequencing_infos_sequencing_info_id) REFERENCES "SequencingInfo" (sequencing_info_id), 
	FOREIGN KEY(panel_info_id) REFERENCES "PanelInfo" (id), 
	FOREIGN KEY(microhaplotypes_detected_id) REFERENCES "MicrohaplotypesDetected" (id), 
	FOREIGN KEY(target_demultiplexed_experiment_samples_id) REFERENCES "DemultiplexedExperimentSamples" (id), 
	FOREIGN KEY(taramp_bioinformatics_infos_tar_amp_bioinformatics_info_id) REFERENCES "TarAmpBioinformaticsInfo" (tar_amp_bioinformatics_info_id), 
	FOREIGN KEY(postprocessing_bioinformatics_infos_id) REFERENCES "BioMethod" (id)
);
CREATE TABLE "RepresentativeMicrohaplotypeSequences" (
	tar_amp_bioinformatics_info_id TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_analysis_name" TEXT, 
	PRIMARY KEY (tar_amp_bioinformatics_info_id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_analysis_name") REFERENCES "PortableMicrohaplotypeObject" (analysis_name)
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
CREATE TABLE "ExperimentInfo" (
	sequencing_info_id TEXT NOT NULL, 
	plate_name TEXT, 
	plate_row TEXT, 
	plate_col INTEGER, 
	specimen_id TEXT NOT NULL, 
	panel_id TEXT NOT NULL, 
	experiment_sample_id TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_analysis_name" TEXT, 
	PRIMARY KEY (experiment_sample_id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_analysis_name") REFERENCES "PortableMicrohaplotypeObject" (analysis_name)
);
CREATE TABLE "SpecimenInfo" (
	specimen_id TEXT NOT NULL, 
	plate_name TEXT, 
	plate_row TEXT, 
	plate_col INTEGER, 
	samp_taxon_id INTEGER NOT NULL, 
	individual_id TEXT, 
	host_taxon_id INTEGER, 
	parasite_density INTEGER, 
	collection_date TEXT NOT NULL, 
	collection_country TEXT NOT NULL, 
	geo_admin1 TEXT, 
	geo_admin2 TEXT, 
	geo_admin3 TEXT, 
	lat_lon TEXT, 
	collector TEXT NOT NULL, 
	samp_store_loc TEXT NOT NULL, 
	samp_collect_device TEXT NOT NULL, 
	project_name TEXT NOT NULL, 
	accession TEXT, 
	sample_comments TEXT, 
	"PortableMicrohaplotypeObject_analysis_name" TEXT, 
	PRIMARY KEY (specimen_id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_analysis_name") REFERENCES "PortableMicrohaplotypeObject" (analysis_name)
);
CREATE TABLE "RepresentativeMicrohaplotypesForTarget" (
	target_id TEXT NOT NULL, 
	"RepresentativeMicrohaplotypeSequences_tar_amp_bioinformatics_info_id" TEXT, 
	PRIMARY KEY (target_id), 
	FOREIGN KEY("RepresentativeMicrohaplotypeSequences_tar_amp_bioinformatics_info_id") REFERENCES "RepresentativeMicrohaplotypeSequences" (tar_amp_bioinformatics_info_id)
);
CREATE TABLE "SpecimenInfo_alternate_identifiers" (
	"SpecimenInfo_specimen_id" TEXT, 
	alternate_identifiers TEXT, 
	PRIMARY KEY ("SpecimenInfo_specimen_id", alternate_identifiers), 
	FOREIGN KEY("SpecimenInfo_specimen_id") REFERENCES "SpecimenInfo" (specimen_id)
);
CREATE TABLE "RepresentativeMicrohaplotypesForTarget_seqs" (
	"RepresentativeMicrohaplotypesForTarget_target_id" TEXT, 
	seqs_id INTEGER NOT NULL, 
	PRIMARY KEY ("RepresentativeMicrohaplotypesForTarget_target_id", seqs_id), 
	FOREIGN KEY("RepresentativeMicrohaplotypesForTarget_target_id") REFERENCES "RepresentativeMicrohaplotypesForTarget" (target_id), 
	FOREIGN KEY(seqs_id) REFERENCES "RepresentativeMicrohaplotypeSequence" (id)
);