-- # Class: "TargetInfo" Description: "Information about a specific target within a genome"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
--     * Slot: gene_id Description: an identifier of the gene, if any, is being covered with this targeted
--     * Slot: PanelInfo_id Description: Autocreated FK slot
--     * Slot: insert_location_id Description: the intended genomic location of the insert of the amplicon (the location between the end of the forward primer and the beginning of the reverse primer)
--     * Slot: forward_primers_id Description: A holder of forward primers associated with this target
--     * Slot: reverse_primers_id Description: A holder of reverse primers associated with this target
-- # Class: "PanelInfo" Description: "information on a panel of targeted amplicon primer pairs"
--     * Slot: id Description: 
--     * Slot: panel_id Description: name of the panel
--     * Slot: target_genome_id Description: the info on the target reference genome for this panel
-- # Class: "RepresentativeHaplotypeSequence" Description: "the representative sequence for a haplotype, similar to a fast(a/q) format"
--     * Slot: id Description: 
--     * Slot: seq Description: the DNA sequence
--     * Slot: haplotype_id Description: name of the haplotype, should be unique to this haplotype
--     * Slot: quality Description: the ansi fastq per base quality score for this sequence, this is optional
--     * Slot: RepresentativeHaplotypeSequences_id Description: Autocreated FK slot
-- # Class: "RepresentativeHaplotypeSequences" Description: "a list of the representative sequence for a haplotypes, similar to a fast(a/q) format"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
-- # Class: "HaplotypesDetected" Description: "the haplotypes detected in a targeted amplicon analysis"
--     * Slot: id Description: 
--     * Slot: sequencing_id Description: the name of the sequencing/wet lab processing steps associated with this run
--     * Slot: bioinformatics_id Description: the name of the bioinformatics processing steps associated with this run
-- # Class: "DemultiplexedSamples" Description: "a list of raw reads counts for each sample for all targets within panel "
--     * Slot: id Description: 
-- # Class: "DemultiplexedTargetsForSample" Description: "a list of raw reads for a sample for all targets within panel"
--     * Slot: sample_id Description: the name of the sample
-- # Class: "DemultiplexedTargetForSample" Description: "the raw read count for a sample for a target"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
--     * Slot: raw_read_count Description: the raw read counts extracted for a target for a sample
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
-- # Class: "PrimerInfo" Description: "information on a primer sequence"
--     * Slot: id Description: 
--     * Slot: seq Description: the DNA sequence
--     * Slot: Primers_id Description: Autocreated FK slot
--     * Slot: location_id Description: what the intended genomic location of the primer is    
-- # Class: "Primers" Description: "A holder of primer sequences"
--     * Slot: id Description: 
-- # Class: "HaplotypesForSample" Description: "Haplotypes detected for a sample for all targets"
--     * Slot: sample_id Description: the name of the sample
-- # Class: "HaplotypeForTarget" Description: "Haplotype detected for a specific target"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
--     * Slot: haplotype_id Description: name of the haplotype, should be unique to this haplotype
--     * Slot: read_count Description: the read count associated with this haplotype
--     * Slot: umi_count Description: the unique molecular identifier (umi) count associated with this haplotype
--     * Slot: HaplotypesForTarget_id Description: Autocreated FK slot
-- # Class: "HaplotypesForTarget" Description: "Haplotypes detected for a specific target"
--     * Slot: id Description: 
--     * Slot: target_id Description: name of the target
-- # Class: "TarAmpBioinformaticsInfo" Description: "the targeted amplicon bioinformatics pipeline"
--     * Slot: tar_amp_bioinformatics_info_id Description: a unique identifier for this bioinformatics info
--     * Slot: demultiplexing_method_id Description: the demultiplexing method used to separate raw reads from barcodes and primer targets
--     * Slot: denoising_method_id Description: the method used to de-noise and/or cluster the raw reads
--     * Slot: population_clustering_method_id Description: the method used to compare clustered/de-noised reads across samples for a target
-- # Class: "BioMethod" Description: "methodology description of a portion of a bioinformatics pipeline"
--     * Slot: id Description: 
--     * Slot: program Description: name of the program used for this portion of the pipeline
--     * Slot: purpose Description: the propose for this method
--     * Slot: program_version Description: a versioning info for the program
--     * Slot: TarAmpBioinformaticsInfo_tar_amp_bioinformatics_info_id Description: Autocreated FK slot
-- # Class: "SequencingInfo" Description: "Information on sequencing info"
--     * Slot: id Description: 
--     * Slot: sequencing_info_id Description: a unique identifier for this sequencing info
--     * Slot: seq_instrument Description: the sequencing instrument used to sequence the run, e.g. ILLUMINA, Illumina MiSeq
--     * Slot: seq_date Description: the date of sequencing, should be YYYY-MM or YYYY-MM-DD
--     * Slot: nucl_acid_ext Description: Link to a reference or kit that describes the recovery of nucleic acids from the sample
--     * Slot: nucl_acid_amp Description: Link to a reference or kit that describes the enzymatic amplification of nucleic acids,
--     * Slot: nucl_acid_date Description: the date of the extraction/amplification
--     * Slot: pcr_cond Description: the method/conditions for PCR, List PCR cycles used to amplify the target
--     * Slot: lib_screen Description: Describe enrichment, screening, or normalization methods applied during amplification or library preparation, e.g. size selection 390bp, diluted to 1 ng DNA/sample
--     * Slot: lib_layout Description: Specify the configuration of reads, e.g. paired-end
--     * Slot: lib_kit Description: Name, version, and applicable cell or cycle numbers for the kit used to prepare libraries and load cells or chips for sequencing. If possible, include a part number, e.g. MiSeq Reagent Kit v3 (150-cycle), MS-102-3001
--     * Slot: seq_center Description: Name of facility where sequencing was performed (lab, core facility, or company)
-- # Class: "SpecimenInfo" Description: "Information on specimen info"
--     * Slot: sample_id Description: the name of the sample
--     * Slot: samp_taxon_id Description: the NCBI taxonomy number of the organism of interest
--     * Slot: host_taxon_id Description: optional the NCBI taxonomy number of the host of the organism
--     * Slot: parasite_density Description: the parasite density in microliters
--     * Slot: collection_date Description: the date of the sample collection
--     * Slot: lat_lon Description: the latitude and longitude of the collection site of the specimen
--     * Slot: collector Description: the name of the primary person managing the specimen
--     * Slot: samp_store_loc Description: the sample store site, address or facility name
--     * Slot: samp_collect_device Description: the way the sample was collected, e.g. whole blood, dried blood spot, etc
--     * Slot: project_name Description: a name of the project under which the sample is organized
--     * Slot: accession Description: ERA/SRA accession number for the sample if it was submitted
--     * Slot: sample_comments Description: any additional comments about the sample
--     * Slot: plate_name Description: a name of plate the specimen was in
--     * Slot: plate_row Description: the row the specimen was in
--     * Slot: plate_col Description: the column the specimen was in 
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "PortableMicrohaplotypeObject" Description: "Information on final results from a targeted amplicon analysis"
--     * Slot: id Description: 
--     * Slot: bioinformatics_info Description: the bioinformatics pipeline info for this project
--     * Slot: sequencing_info_id Description: the sequencing info for this project
--     * Slot: panel_info_id Description: the info on the panel used for this project
--     * Slot: haplotypes_detected_id Description: the haplotypes detected in this projects
--     * Slot: target_demultiplexed_samples_id Description: the raw demultiplex target counts for each sample  
-- # Class: "RepresentativeHaplotypeSequence_alt_annotations" Description: ""
--     * Slot: RepresentativeHaplotypeSequence_id Description: Autocreated FK slot
--     * Slot: alt_annotations Description: a list of additional annotations associated with this haplotype, e.g. wildtype, amino acid changes etc
-- # Class: "HaplotypesDetected_samples" Description: ""
--     * Slot: HaplotypesDetected_id Description: Autocreated FK slot
--     * Slot: samples_sample_id Description: a list of the haplotypes detected for a sample for various targets 
-- # Class: "DemultiplexedSamples_demultiplexed_samples" Description: ""
--     * Slot: DemultiplexedSamples_id Description: Autocreated FK slot
--     * Slot: demultiplexed_samples_sample_id Description: a list of the samples with the number of raw reads extracted 
-- # Class: "DemultiplexedTargetsForSample_demultiplexed_targets" Description: ""
--     * Slot: DemultiplexedTargetsForSample_sample_id Description: Autocreated FK slot
--     * Slot: demultiplexed_targets_id Description: a list of the targets extracted for a sample 
-- # Class: "HaplotypesForSample_target_results" Description: ""
--     * Slot: HaplotypesForSample_sample_id Description: Autocreated FK slot
--     * Slot: target_results_id Description: a list of the haplotypes detected for a list of targets
-- # Class: "BioMethod_additional_argument" Description: ""
--     * Slot: BioMethod_id Description: Autocreated FK slot
--     * Slot: additional_argument Description: any additional arguments that differ from the default
-- # Class: "SpecimenInfo_alternate_identifiers" Description: ""
--     * Slot: SpecimenInfo_sample_id Description: Autocreated FK slot
--     * Slot: alternate_identifiers Description: a list of optional alternative names for the samples
-- # Class: "PortableMicrohaplotypeObject_representative_haplotype_sequences" Description: ""
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
--     * Slot: representative_haplotype_sequences_id Description: a list of the representative sequences for the results for this project

CREATE TABLE "RepresentativeHaplotypeSequences" (
	id INTEGER NOT NULL, 
	target_id TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "HaplotypesDetected" (
	id INTEGER NOT NULL, 
	sequencing_id TEXT NOT NULL, 
	bioinformatics_id TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "DemultiplexedSamples" (
	id INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "DemultiplexedTargetsForSample" (
	sample_id TEXT NOT NULL, 
	PRIMARY KEY (sample_id)
);
CREATE TABLE "DemultiplexedTargetForSample" (
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
	PRIMARY KEY (id)
);
CREATE TABLE "Primers" (
	id INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "HaplotypesForSample" (
	sample_id TEXT NOT NULL, 
	PRIMARY KEY (sample_id)
);
CREATE TABLE "HaplotypesForTarget" (
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
	id INTEGER NOT NULL, 
	sequencing_info_id TEXT NOT NULL, 
	seq_instrument TEXT NOT NULL, 
	seq_date TEXT NOT NULL, 
	nucl_acid_ext TEXT NOT NULL, 
	nucl_acid_amp TEXT NOT NULL, 
	nucl_acid_date TEXT NOT NULL, 
	pcr_cond TEXT NOT NULL, 
	lib_screen TEXT NOT NULL, 
	lib_layout TEXT NOT NULL, 
	lib_kit TEXT NOT NULL, 
	seq_center TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "PanelInfo" (
	id INTEGER NOT NULL, 
	panel_id TEXT NOT NULL, 
	target_genome_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(target_genome_id) REFERENCES "GenomeInfo" (id)
);
CREATE TABLE "RepresentativeHaplotypeSequence" (
	id INTEGER NOT NULL, 
	seq TEXT NOT NULL, 
	haplotype_id TEXT NOT NULL, 
	quality TEXT, 
	"RepresentativeHaplotypeSequences_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("RepresentativeHaplotypeSequences_id") REFERENCES "RepresentativeHaplotypeSequences" (id)
);
CREATE TABLE "PrimerInfo" (
	id INTEGER NOT NULL, 
	seq TEXT NOT NULL, 
	"Primers_id" INTEGER, 
	location_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY("Primers_id") REFERENCES "Primers" (id), 
	FOREIGN KEY(location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "HaplotypeForTarget" (
	id INTEGER NOT NULL, 
	target_id TEXT NOT NULL, 
	haplotype_id TEXT NOT NULL, 
	read_count FLOAT NOT NULL, 
	umi_count FLOAT, 
	"HaplotypesForTarget_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("HaplotypesForTarget_id") REFERENCES "HaplotypesForTarget" (id)
);
CREATE TABLE "HaplotypesDetected_samples" (
	"HaplotypesDetected_id" INTEGER, 
	samples_sample_id TEXT NOT NULL, 
	PRIMARY KEY ("HaplotypesDetected_id", samples_sample_id), 
	FOREIGN KEY("HaplotypesDetected_id") REFERENCES "HaplotypesDetected" (id), 
	FOREIGN KEY(samples_sample_id) REFERENCES "HaplotypesForSample" (sample_id)
);
CREATE TABLE "DemultiplexedSamples_demultiplexed_samples" (
	"DemultiplexedSamples_id" INTEGER, 
	demultiplexed_samples_sample_id TEXT NOT NULL, 
	PRIMARY KEY ("DemultiplexedSamples_id", demultiplexed_samples_sample_id), 
	FOREIGN KEY("DemultiplexedSamples_id") REFERENCES "DemultiplexedSamples" (id), 
	FOREIGN KEY(demultiplexed_samples_sample_id) REFERENCES "DemultiplexedTargetsForSample" (sample_id)
);
CREATE TABLE "DemultiplexedTargetsForSample_demultiplexed_targets" (
	"DemultiplexedTargetsForSample_sample_id" TEXT, 
	demultiplexed_targets_id INTEGER NOT NULL, 
	PRIMARY KEY ("DemultiplexedTargetsForSample_sample_id", demultiplexed_targets_id), 
	FOREIGN KEY("DemultiplexedTargetsForSample_sample_id") REFERENCES "DemultiplexedTargetsForSample" (sample_id), 
	FOREIGN KEY(demultiplexed_targets_id) REFERENCES "DemultiplexedTargetForSample" (id)
);
CREATE TABLE "HaplotypesForSample_target_results" (
	"HaplotypesForSample_sample_id" TEXT, 
	target_results_id INTEGER NOT NULL, 
	PRIMARY KEY ("HaplotypesForSample_sample_id", target_results_id), 
	FOREIGN KEY("HaplotypesForSample_sample_id") REFERENCES "HaplotypesForSample" (sample_id), 
	FOREIGN KEY(target_results_id) REFERENCES "HaplotypesForTarget" (id)
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
	insert_location_id INTEGER NOT NULL, 
	forward_primers_id INTEGER NOT NULL, 
	reverse_primers_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PanelInfo_id") REFERENCES "PanelInfo" (id), 
	FOREIGN KEY(insert_location_id) REFERENCES "GenomicLocation" (id), 
	FOREIGN KEY(forward_primers_id) REFERENCES "Primers" (id), 
	FOREIGN KEY(reverse_primers_id) REFERENCES "Primers" (id)
);
CREATE TABLE "PortableMicrohaplotypeObject" (
	id INTEGER NOT NULL, 
	bioinformatics_info TEXT NOT NULL, 
	sequencing_info_id INTEGER NOT NULL, 
	panel_info_id INTEGER NOT NULL, 
	haplotypes_detected_id INTEGER NOT NULL, 
	target_demultiplexed_samples_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(bioinformatics_info) REFERENCES "TarAmpBioinformaticsInfo" (tar_amp_bioinformatics_info_id), 
	FOREIGN KEY(sequencing_info_id) REFERENCES "SequencingInfo" (id), 
	FOREIGN KEY(panel_info_id) REFERENCES "PanelInfo" (id), 
	FOREIGN KEY(haplotypes_detected_id) REFERENCES "HaplotypesDetected" (id), 
	FOREIGN KEY(target_demultiplexed_samples_id) REFERENCES "DemultiplexedSamples" (id)
);
CREATE TABLE "RepresentativeHaplotypeSequence_alt_annotations" (
	"RepresentativeHaplotypeSequence_id" INTEGER, 
	alt_annotations TEXT, 
	PRIMARY KEY ("RepresentativeHaplotypeSequence_id", alt_annotations), 
	FOREIGN KEY("RepresentativeHaplotypeSequence_id") REFERENCES "RepresentativeHaplotypeSequence" (id)
);
CREATE TABLE "SpecimenInfo" (
	sample_id TEXT NOT NULL, 
	samp_taxon_id INTEGER NOT NULL, 
	host_taxon_id INTEGER, 
	parasite_density INTEGER, 
	collection_date TEXT NOT NULL, 
	lat_lon TEXT NOT NULL, 
	collector TEXT NOT NULL, 
	samp_store_loc TEXT NOT NULL, 
	samp_collect_device TEXT NOT NULL, 
	project_name TEXT NOT NULL, 
	accession TEXT, 
	sample_comments TEXT, 
	plate_name TEXT NOT NULL, 
	plate_row TEXT NOT NULL, 
	plate_col TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (sample_id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "PortableMicrohaplotypeObject_representative_haplotype_sequences" (
	"PortableMicrohaplotypeObject_id" INTEGER, 
	representative_haplotype_sequences_id INTEGER NOT NULL, 
	PRIMARY KEY ("PortableMicrohaplotypeObject_id", representative_haplotype_sequences_id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id), 
	FOREIGN KEY(representative_haplotype_sequences_id) REFERENCES "RepresentativeHaplotypeSequences" (id)
);
CREATE TABLE "SpecimenInfo_alternate_identifiers" (
	"SpecimenInfo_sample_id" TEXT, 
	alternate_identifiers TEXT, 
	PRIMARY KEY ("SpecimenInfo_sample_id", alternate_identifiers), 
	FOREIGN KEY("SpecimenInfo_sample_id") REFERENCES "SpecimenInfo" (sample_id)
);