-- # Class: "MarkerOfInterest" Description: "A specific genomic location of interest, e.g. drug resistance, or other phenotypical marker"
--     * Slot: id Description: 
--     * Slot: marker_location_id Description: the genomic location
-- # Class: "TargetInfo" Description: "Information about a specific target within a genome"
--     * Slot: id Description: 
--     * Slot: gene_name Description: an identifier of the gene, if any, is being covered with this targeted
--     * Slot: target_name Description: an identifier for this target
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
--     * Slot: insert_location_id Description: the intended genomic location of the insert of the amplicon (the location between the end of the forward primer and the beginning of the reverse primer)
--     * Slot: forward_primer_id Description: the forward primer associated with this target
--     * Slot: reverse_primer_id Description: the reverse primer associated with this target
-- # Class: "ReactionInfo" Description: "information on a panel of targeted amplicon primer pairs"
--     * Slot: id Description: 
--     * Slot: reaction_name Description: a name for this reaction
-- # Class: "PanelInfo" Description: "information on a panel of targeted amplicon primer pairs"
--     * Slot: id Description: 
--     * Slot: panel_name Description: a name for the panel
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "Pseudocigar" Description: "information on pseudocigar for a sequence"
--     * Slot: id Description: 
--     * Slot: pseudocigar_seq Description: the pseudocigar itself
--     * Slot: pseudocigar_generation_description Description: a description of how the pseudocigar information was generated 
--     * Slot: ref_loc_id Description: the genomic location the pseudocigar is in reference to
-- # Class: "RepresentativeMicrohaplotype" Description: "the representative sequence for a microhaplotype, similar to a fast(a/q) format"
--     * Slot: id Description: 
--     * Slot: seq Description: the DNA sequence
--     * Slot: microhaplotype_name Description: an optional name for this microhaplotype
--     * Slot: quality Description: the ansi fastq per base quality score for this sequence, this is optional
--     * Slot: pseudocigar_id Description: the pseudocigar of the haplotype
-- # Class: "MaskingInfo" Description: "information about a subsegment of the sequence that should be masked"
--     * Slot: id Description: 
--     * Slot: seq_start Description: the start of the masking
--     * Slot: seq_segment_size Description: the size of the masking
--     * Slot: replacement_size Description: the size of replacement mask
--     * Slot: masking_generation_description Description: a description of how the masking information was generated 
-- # Class: "RepresentativeMicrohaplotypes" Description: "a collection of representative sequences for microhaplotypes for all targets"
--     * Slot: id Description: 
-- # Class: "RepresentativeMicrohaplotypesForTarget" Description: "a list of the representative sequence for a microhaplotypes, similar to a fast(a/q) format"
--     * Slot: target_id Description: the index into the target_info list
--     * Slot: RepresentativeMicrohaplotypes_id Description: Autocreated FK slot
--     * Slot: mhap_location_id Description: a genomic location that was analyzed for this target info, this allows listing location that may be different from the full target location (e.g 1 base in from the full) 
-- # Class: "DetectedMicrohaplotypes" Description: "the microhaplotypes detected in a targeted amplicon analysis"
--     * Slot: id Description: 
--     * Slot: bioinformatics_run_id Description: the index into bioinformatics_run_info list
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "GenomeInfo" Description: "information on a genome"
--     * Slot: id Description: 
--     * Slot: name Description: name of the genome
--     * Slot: genome_version Description: the genome version
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
--     * Slot: alt_seq Description: a possible alternative sequence of this genomic location
--     * Slot: RepresentativeMicrohaplotype_id Description: Autocreated FK slot
-- # Class: "ProteinVariant" Description: "information on a variant in protein sequence"
--     * Slot: id Description: 
--     * Slot: gene_name Description: an identifier of the gene, if any, is being covered with this targeted
--     * Slot: alternative_gene_name Description: an alternative gene name
--     * Slot: RepresentativeMicrohaplotype_id Description: Autocreated FK slot
--     * Slot: protein_location_id Description: the position within the protein, the chromosome in this case would be the transcript name
--     * Slot: codon_genomic_location_id Description: the position within the genomic sequence of the codon
-- # Class: "PrimerInfo" Description: "information on a primer sequence"
--     * Slot: id Description: 
--     * Slot: seq Description: the DNA sequence
--     * Slot: location_id Description: what the intended genomic location of the primer is
-- # Class: "DetectedMicrohaplotypesForSample" Description: "Microhaplotypes detected for a sample for all targets"
--     * Slot: id Description: 
--     * Slot: library_sample_id Description: the index into the library_sample_info list
--     * Slot: DetectedMicrohaplotypes_id Description: Autocreated FK slot
-- # Class: "MicrohaplotypeForTarget" Description: "Microhaplotype detected for a specific target"
--     * Slot: id Description: 
--     * Slot: mhap_id Description: the index for a microhaplotype for a target in the representative_microhaplotypes list, e.g. representative_microhaplotypes[mhaps_target_id][mhap_id]
--     * Slot: reads Description: the read count associated with this microhaplotype
--     * Slot: umis Description: the unique molecular identifier (umi) count associated with this microhaplotype
--     * Slot: DetectedMicrohaplotypesForTarget_id Description: Autocreated FK slot
-- # Class: "DetectedMicrohaplotypesForTarget" Description: "Microhaplotypes detected for a specific target"
--     * Slot: id Description: 
--     * Slot: mhaps_target_id Description: the index for a target in the representative_microhaplotypes list
--     * Slot: DetectedMicrohaplotypesForSample_id Description: Autocreated FK slot
-- # Class: "BioinformaticsMethodInfo" Description: "the targeted amplicon bioinformatics methods used to generate the amplicon data in this PMO"
--     * Slot: id Description: 
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "BioMethod" Description: "methodology description of a portion of a bioinformatics pipeline"
--     * Slot: id Description: 
--     * Slot: program_version Description: the version of generation method, should be in the format of v[MAJOR].[MINOR].[PATCH]
--     * Slot: program Description: name of the program used for this portion of the pipeline
--     * Slot: program_description Description: a short description of what this method does
--     * Slot: program_url Description: a url pointing to code base of a program, e.g. a github link
--     * Slot: BioinformaticsMethodInfo_id Description: Autocreated FK slot
-- # Class: "PlateInfo" Description: "Information about a plate location in a standard 96 well plate"
--     * Slot: id Description: 
--     * Slot: plate_name Description: a name of plate the specimen was in
--     * Slot: plate_row Description: the row the specimen was in
--     * Slot: plate_col Description: the column the specimen was in
-- # Class: "LibrarySampleInfo" Description: "Information about a specific amplification and sequencing of a specimen"
--     * Slot: sequencing_info_id Description: the index into the sequencing_info list
--     * Slot: specimen_id Description: the index into the specimen_info list
--     * Slot: panel_id Description: the index into the panel_info list
--     * Slot: fastqs_loc Description: the location (url or filename path) of the fastqs for a library run
--     * Slot: run_accession Description: ERA/SRA run accession number for the sample if it was submitted
--     * Slot: experiment_accession Description: ERA/SRA experiment accession number for the sample if it was submitted
--     * Slot: library_sample_name Description: a unique identifier for this sequence/amplification run on a specimen_name
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
--     * Slot: library_prep_plate_info_id Description: plate location of where library was prepared for sequencing 
-- # Class: "SequencingInfo" Description: "Information on sequencing info"
--     * Slot: id Description: 
--     * Slot: sequencing_info_name Description: a name of for the sequencing done, e.g. batch1
--     * Slot: seq_platform Description: the sequencing technology used to sequence the run, e.g. ILLUMINA, NANOPORE, PACBIO
--     * Slot: seq_instrument_model Description: the sequencing instrument model used to sequence the run, e.g. NextSeq 2000, MinION, Revio
--     * Slot: seq_date Description: the date of sequencing, should be YYYY-MM or YYYY-MM-DD
--     * Slot: nucl_acid_ext Description: Link to a reference or kit that describes the recovery of nucleic acids from the sample
--     * Slot: nucl_acid_amp Description: Link to a reference or kit that describes the enzymatic amplification of nucleic acids
--     * Slot: nucl_acid_ext_date Description: the date of the nucleoacid extraction
--     * Slot: nucl_acid_amp_date Description: the date of the nucleoacid amplification
--     * Slot: pcr_cond Description: the method/conditions for PCR, List PCR cycles used to amplify the target
--     * Slot: library_screen Description: Describe enrichment, screening, or normalization methods applied during amplification or library preparation, e.g. size selection 390bp, diluted to 1 ng DNA/sample
--     * Slot: library_kit Description: Name, version, and applicable cell or cycle numbers for the kit used to prepare libraries and load cells or chips for sequencing. If possible, include a part number, e.g. MiSeq Reagent Kit v3 (150-cycle), MS-102-3001
--     * Slot: library_layout Description: Specify the configuration of reads, e.g. paired-end, single
--     * Slot: library_strategy Description: what the nuceloacid sequencing/amplification strategy was (common names are AMPLICON, WGS)
--     * Slot: library_source Description: Source of amplification material (common names GENOMIC, TRANSCRIPTOMIC)
--     * Slot: library_selection Description: how amplification was done (common are PCR=Source material was selected by designed primers, RANDOM =Random selection by shearing or other method)
--     * Slot: seq_center Description: Name of facility where sequencing was performed (lab, core facility, or company)
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "ParasiteDensity" Description: "method and value of determined parasite density"
--     * Slot: id Description: 
--     * Slot: parasite_density_method Description: the method of how this density was obtained
--     * Slot: parasite_density Description: the density in microliters
--     * Slot: date_measured Description: the date the qpcr was performed, can be YYYY, YYYY-MM, or YYYY-MM-DD
--     * Slot: density_method_comments Description: additional comments about how the density was performed
--     * Slot: LibrarySampleInfo_library_sample_name Description: Autocreated FK slot
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
-- # Class: "ProjectInfo" Description: "Information on project info"
--     * Slot: project_name Description: a name for the project, should be unique if multiple projects listed
--     * Slot: project_description Description: a short description of the project
--     * Slot: project_type Description: the type of project conducted, e.g. TES vs surveillance vs transmission
--     * Slot: project_collector_chief_scientist Description: can be collection of names separated by a semicolon if multiple people involved or can just be the name of the primary person managing the specimen
--     * Slot: BioProject_accession Description: an SRA bioproject accession e.g. PRJNA33823
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "TravelInfo" Description: "Information on travel info"
--     * Slot: id Description: 
--     * Slot: lat_lon Description: the latitude and longitude of a specific site
--     * Slot: geo_admin1 Description: geographical admin level 1, the secondary large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin2 Description: geographical admin level 2, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin3 Description: geographical admin level 3, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: travel_country Description: the name of country, would be the same as admin level 0
--     * Slot: travel_start_date Description: the date of the start of travel, can be approximate, should be YYYY-MM or YYYY-MM-DD (preferred)
--     * Slot: travel_end_date Description: the date of the end of travel, can be approximate, should be YYYY-MM or YYYY-MM-DD (preferred)
--     * Slot: bed_net_usage Description: approximate usage of bed net while traveling, 1 = 100% nights with bed net, 0 = 0% no bed net usage
-- # Class: "SpecimenInfo" Description: "Information on specimen info"
--     * Slot: lat_lon Description: the latitude and longitude of a specific site
--     * Slot: geo_admin1 Description: geographical admin level 1, the secondary large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin2 Description: geographical admin level 2, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: geo_admin3 Description: geographical admin level 3, the third large demarcation of a nation (nation = admin level 0)
--     * Slot: specimen_name Description: an identifier for the specimen, should be unique within this sample set
--     * Slot: host_subject_id Description: an identifier for the individual a specimen was collected from
--     * Slot: host_taxon_id Description: the NCBI taxonomy number of the host that the specimen was collected from
--     * Slot: host_sex Description: if specimen is from a person, the sex listed for that person
--     * Slot: specimen_accession Description: if specimen is deposited in a database, what accession is it associated with
--     * Slot: gravid Description: whether host specimen is currently pregnant
--     * Slot: blood_meal Description: whether host specimen has had a recent blood meal  
--     * Slot: gravidity Description: the gravidity of the specimen host (number of previous pregnancies)
--     * Slot: collection_date Description: the date of the specimen collection, can be YYYY, YYYY-MM, or YYYY-MM-DD
--     * Slot: host_age Description: if specimen is from a person, the age in years of the person, can be float value so for 3 month old put 0.25
--     * Slot: collection_country Description: the name of country collected in, would be the same as admin level 0
--     * Slot: specimen_store_loc Description: the specimen store site, address or facility name
--     * Slot: specimen_collect_device Description: the way the specimen was collected, e.g. whole blood, dried blood spot
--     * Slot: specimen_type Description: what type of specimen this is, e.g. negative_control, positive_control, field_sample
--     * Slot: project_id Description: the index into the project_info list
--     * Slot: has_travel_out_six_month Description: has travelled out from local region in the last six months
--     * Slot: env_medium Description: the environment medium from which the specimen was collected from
--     * Slot: env_local_scale Description: the local environment from which the specimen was collected, e.g. jungle, urban, rural
--     * Slot: env_broad_scale Description: the broad environment from which the specimen was collected, e.g. highlands, lowlands, mountainous region
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
--     * Slot: storage_plate_info_id Description: plate location of where specimen is stored if stored in a plate 
-- # Class: "BioinformaticsRunInfo" Description: "Information about the pipeline run that generated some of the microhaplotype detected and reads_by_stage"
--     * Slot: id Description: 
--     * Slot: bioinformatics_methods_id Description: the index into the bioinformatics_methods_info list
--     * Slot: run_date Description: the date when the run was done, should be YYYY-MM-DD
--     * Slot: bioinformatics_run_name Description: a name to for this run, needs to be unique to each run 
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
--     * Slot: reads Description: the read counts for this stage
--     * Slot: stage Description: the stage of the pipeline, e.g. demultiplexed, denoised, etc
--     * Slot: ReadCountsByStageForTarget_id Description: Autocreated FK slot
-- # Class: "ReadCountsByStageForTarget" Description: "Information on the reads counts at several stages of a pipeline for a target"
--     * Slot: id Description: 
--     * Slot: target_id Description: the index into the target_info list
--     * Slot: ReadCountsByStageForLibrarySample_id Description: Autocreated FK slot
-- # Class: "ReadCountsByStageForLibrarySample" Description: "Information on the reads counts at several stages of a pipeline for a library_sample"
--     * Slot: id Description: 
--     * Slot: library_sample_id Description: the index into the library_sample_info list
--     * Slot: total_raw_count Description: the raw counts off the sequencing machine that a sample began with
--     * Slot: ReadCountsByStage_id Description: Autocreated FK slot
-- # Class: "ReadCountsByStage" Description: "Information on the reads counts at several stages of a pipeline"
--     * Slot: id Description: 
--     * Slot: bioinformatics_run_id Description: the index into bioinformatics_run_info list
--     * Slot: PortableMicrohaplotypeObject_id Description: Autocreated FK slot
-- # Class: "PortableMicrohaplotypeObject" Description: "Information on final results from a targeted amplicon analysis"
--     * Slot: id Description: 
--     * Slot: representative_microhaplotypes_id Description: a list of the information on the representative microhaplotypes
--     * Slot: pmo_header_id Description: the PMO information for this file including version etc
-- # Class: "MarkerOfInterest_associations" Description: ""
--     * Slot: MarkerOfInterest_id Description: Autocreated FK slot
--     * Slot: associations Description: a list of associations with this marker, e.g. SP resistance, etc
-- # Class: "TargetInfo_markers_of_interest" Description: ""
--     * Slot: TargetInfo_id Description: Autocreated FK slot
--     * Slot: markers_of_interest_id Description: a list of covered markers of interest
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
--     * Slot: alt_annotations Description: a list of additional annotations associated with this microhaplotype, e.g. wildtype
-- # Class: "RepresentativeMicrohaplotypesForTarget_microhaplotypes" Description: ""
--     * Slot: RepresentativeMicrohaplotypesForTarget_target_id Description: Autocreated FK slot
--     * Slot: microhaplotypes_id Description: a list of the microhaplotypes detected for a target
-- # Class: "GenomeInfo_taxon_id" Description: ""
--     * Slot: GenomeInfo_id Description: Autocreated FK slot
--     * Slot: taxon_id Description: the NCBI taxonomy number, can be a list of values
-- # Class: "GenomeInfo_chromosomes" Description: ""
--     * Slot: GenomeInfo_id Description: Autocreated FK slot
--     * Slot: chromosomes Description: a list of chromosomes found within this genome
-- # Class: "BioMethod_additional_argument" Description: ""
--     * Slot: BioMethod_id Description: Autocreated FK slot
--     * Slot: additional_argument Description: any additional arguments that differ from the default
-- # Class: "LibrarySampleInfo_alternate_identifiers" Description: ""
--     * Slot: LibrarySampleInfo_library_sample_name Description: Autocreated FK slot
--     * Slot: alternate_identifiers Description: a list of optional alternative names
-- # Class: "ProjectInfo_project_contributors" Description: ""
--     * Slot: ProjectInfo_project_name Description: Autocreated FK slot
--     * Slot: project_contributors Description: a list of collaborators who contributed to this project
-- # Class: "SpecimenInfo_alternate_identifiers" Description: ""
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
--     * Slot: alternate_identifiers Description: a list of optional alternative names
-- # Class: "SpecimenInfo_specimen_taxon_id" Description: ""
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
--     * Slot: specimen_taxon_id Description: the NCBI taxonomy number of the organism in specimen, can list multiple if a mixed sample
-- # Class: "SpecimenInfo_specimen_comments" Description: ""
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
--     * Slot: specimen_comments Description: any additional comments about the specimen
-- # Class: "SpecimenInfo_travel_out_six_month" Description: ""
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
--     * Slot: travel_out_six_month_id Description: Specification of the countries travelled in the last six months; can include multiple travels
-- # Class: "SpecimenInfo_drug_usage" Description: ""
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
--     * Slot: drug_usage Description: Any drug used by subject and the frequency of usage; can include multiple drugs used
-- # Class: "SpecimenInfo_treatment_status" Description: ""
--     * Slot: SpecimenInfo_specimen_name Description: Autocreated FK slot
--     * Slot: treatment_status Description: If person has been treated with drugs, what was the treatment outcome

CREATE TABLE "ReactionInfo" (
	id INTEGER NOT NULL, 
	reaction_name TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "Pseudocigar" (
	id INTEGER NOT NULL, 
	pseudocigar_seq TEXT NOT NULL, 
	pseudocigar_generation_description TEXT, 
	ref_loc_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(ref_loc_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "RepresentativeMicrohaplotype" (
	id INTEGER NOT NULL, 
	seq TEXT NOT NULL, 
	microhaplotype_name TEXT, 
	quality TEXT, 
	pseudocigar_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(pseudocigar_id) REFERENCES "Pseudocigar" (id)
);
CREATE TABLE "MaskingInfo" (
	id INTEGER NOT NULL, 
	seq_start INTEGER NOT NULL, 
	seq_segment_size INTEGER NOT NULL, 
	replacement_size INTEGER NOT NULL, 
	masking_generation_description TEXT, 
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
	alt_seq TEXT, 
	"RepresentativeMicrohaplotype_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("RepresentativeMicrohaplotype_id") REFERENCES "RepresentativeMicrohaplotype" (id)
);
CREATE TABLE "PlateInfo" (
	id INTEGER NOT NULL, 
	plate_name TEXT, 
	plate_row TEXT, 
	plate_col INTEGER, 
	PRIMARY KEY (id)
);
CREATE TABLE "TravelInfo" (
	id INTEGER NOT NULL, 
	lat_lon TEXT, 
	geo_admin1 TEXT, 
	geo_admin2 TEXT, 
	geo_admin3 TEXT, 
	travel_country TEXT NOT NULL, 
	travel_start_date TEXT NOT NULL, 
	travel_end_date TEXT NOT NULL, 
	bed_net_usage FLOAT, 
	PRIMARY KEY (id)
);
CREATE TABLE "PmoGenerationMethod" (
	id INTEGER NOT NULL, 
	program_version TEXT NOT NULL, 
	program_name TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "MarkerOfInterest" (
	id INTEGER NOT NULL, 
	marker_location_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(marker_location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "RepresentativeMicrohaplotypesForTarget" (
	target_id INTEGER NOT NULL, 
	"RepresentativeMicrohaplotypes_id" INTEGER, 
	mhap_location_id INTEGER, 
	PRIMARY KEY (target_id), 
	FOREIGN KEY("RepresentativeMicrohaplotypes_id") REFERENCES "RepresentativeMicrohaplotypes" (id), 
	FOREIGN KEY(mhap_location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "ProteinVariant" (
	id INTEGER NOT NULL, 
	gene_name TEXT, 
	alternative_gene_name TEXT, 
	"RepresentativeMicrohaplotype_id" INTEGER, 
	protein_location_id INTEGER NOT NULL, 
	codon_genomic_location_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("RepresentativeMicrohaplotype_id") REFERENCES "RepresentativeMicrohaplotype" (id), 
	FOREIGN KEY(protein_location_id) REFERENCES "GenomicLocation" (id), 
	FOREIGN KEY(codon_genomic_location_id) REFERENCES "GenomicLocation" (id)
);
CREATE TABLE "PrimerInfo" (
	id INTEGER NOT NULL, 
	seq TEXT NOT NULL, 
	location_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(location_id) REFERENCES "GenomicLocation" (id)
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
CREATE TABLE "PortableMicrohaplotypeObject" (
	id INTEGER NOT NULL, 
	representative_microhaplotypes_id INTEGER NOT NULL, 
	pmo_header_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(representative_microhaplotypes_id) REFERENCES "RepresentativeMicrohaplotypes" (id), 
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
	gene_name TEXT, 
	target_name TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	insert_location_id INTEGER, 
	forward_primer_id INTEGER NOT NULL, 
	reverse_primer_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id), 
	FOREIGN KEY(insert_location_id) REFERENCES "GenomicLocation" (id), 
	FOREIGN KEY(forward_primer_id) REFERENCES "PrimerInfo" (id), 
	FOREIGN KEY(reverse_primer_id) REFERENCES "PrimerInfo" (id)
);
CREATE TABLE "PanelInfo" (
	id INTEGER NOT NULL, 
	panel_name TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "DetectedMicrohaplotypes" (
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
	url TEXT NOT NULL, 
	gff_url TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "BioinformaticsMethodInfo" (
	id INTEGER NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "LibrarySampleInfo" (
	sequencing_info_id INTEGER NOT NULL, 
	specimen_id INTEGER NOT NULL, 
	panel_id INTEGER NOT NULL, 
	fastqs_loc TEXT, 
	run_accession TEXT, 
	experiment_accession TEXT, 
	library_sample_name TEXT NOT NULL, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	library_prep_plate_info_id INTEGER, 
	PRIMARY KEY (library_sample_name), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id), 
	FOREIGN KEY(library_prep_plate_info_id) REFERENCES "PlateInfo" (id)
);
CREATE TABLE "SequencingInfo" (
	id INTEGER NOT NULL, 
	sequencing_info_name TEXT NOT NULL, 
	seq_platform TEXT NOT NULL, 
	seq_instrument_model TEXT NOT NULL, 
	seq_date TEXT, 
	nucl_acid_ext TEXT, 
	nucl_acid_amp TEXT, 
	nucl_acid_ext_date TEXT, 
	nucl_acid_amp_date TEXT, 
	pcr_cond TEXT, 
	library_screen TEXT, 
	library_kit TEXT, 
	library_layout TEXT NOT NULL, 
	library_strategy TEXT NOT NULL, 
	library_source TEXT NOT NULL, 
	library_selection TEXT NOT NULL, 
	seq_center TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "ProjectInfo" (
	project_name TEXT NOT NULL, 
	project_description TEXT NOT NULL, 
	project_type TEXT, 
	project_collector_chief_scientist TEXT, 
	"BioProject_accession" TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	PRIMARY KEY (project_name), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id)
);
CREATE TABLE "SpecimenInfo" (
	lat_lon TEXT, 
	geo_admin1 TEXT, 
	geo_admin2 TEXT, 
	geo_admin3 TEXT, 
	specimen_name TEXT NOT NULL, 
	host_subject_id INTEGER, 
	host_taxon_id INTEGER NOT NULL, 
	host_sex TEXT, 
	specimen_accession TEXT, 
	gravid BOOLEAN, 
	blood_meal BOOLEAN, 
	gravidity INTEGER, 
	collection_date TEXT NOT NULL, 
	host_age FLOAT, 
	collection_country TEXT NOT NULL, 
	specimen_store_loc TEXT, 
	specimen_collect_device TEXT, 
	specimen_type TEXT, 
	project_id INTEGER NOT NULL, 
	has_travel_out_six_month BOOLEAN, 
	env_medium TEXT, 
	env_local_scale TEXT, 
	env_broad_scale TEXT, 
	"PortableMicrohaplotypeObject_id" INTEGER, 
	storage_plate_info_id INTEGER, 
	PRIMARY KEY (specimen_name), 
	FOREIGN KEY("PortableMicrohaplotypeObject_id") REFERENCES "PortableMicrohaplotypeObject" (id), 
	FOREIGN KEY(storage_plate_info_id) REFERENCES "PlateInfo" (id)
);
CREATE TABLE "BioinformaticsRunInfo" (
	id INTEGER NOT NULL, 
	bioinformatics_methods_id INTEGER NOT NULL, 
	run_date TEXT, 
	bioinformatics_run_name TEXT NOT NULL, 
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
CREATE TABLE "DetectedMicrohaplotypesForSample" (
	id INTEGER NOT NULL, 
	library_sample_id INTEGER NOT NULL, 
	"DetectedMicrohaplotypes_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("DetectedMicrohaplotypes_id") REFERENCES "DetectedMicrohaplotypes" (id)
);
CREATE TABLE "BioMethod" (
	id INTEGER NOT NULL, 
	program_version TEXT NOT NULL, 
	program TEXT NOT NULL, 
	program_description TEXT, 
	program_url TEXT, 
	"BioinformaticsMethodInfo_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("BioinformaticsMethodInfo_id") REFERENCES "BioinformaticsMethodInfo" (id)
);
CREATE TABLE "ParasiteDensity" (
	id INTEGER NOT NULL, 
	parasite_density_method TEXT NOT NULL, 
	parasite_density FLOAT NOT NULL, 
	date_measured TEXT, 
	density_method_comments TEXT, 
	"LibrarySampleInfo_library_sample_name" TEXT, 
	"SpecimenInfo_specimen_name" TEXT, 
	PRIMARY KEY (id), 
	FOREIGN KEY("LibrarySampleInfo_library_sample_name") REFERENCES "LibrarySampleInfo" (library_sample_name), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "ReadCountsByStageForLibrarySample" (
	id INTEGER NOT NULL, 
	library_sample_id INTEGER NOT NULL, 
	total_raw_count INTEGER NOT NULL, 
	"ReadCountsByStage_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("ReadCountsByStage_id") REFERENCES "ReadCountsByStage" (id)
);
CREATE TABLE "TargetInfo_markers_of_interest" (
	"TargetInfo_id" INTEGER, 
	markers_of_interest_id INTEGER, 
	PRIMARY KEY ("TargetInfo_id", markers_of_interest_id), 
	FOREIGN KEY("TargetInfo_id") REFERENCES "TargetInfo" (id), 
	FOREIGN KEY(markers_of_interest_id) REFERENCES "MarkerOfInterest" (id)
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
CREATE TABLE "GenomeInfo_taxon_id" (
	"GenomeInfo_id" INTEGER, 
	taxon_id INTEGER NOT NULL, 
	PRIMARY KEY ("GenomeInfo_id", taxon_id), 
	FOREIGN KEY("GenomeInfo_id") REFERENCES "GenomeInfo" (id)
);
CREATE TABLE "GenomeInfo_chromosomes" (
	"GenomeInfo_id" INTEGER, 
	chromosomes TEXT, 
	PRIMARY KEY ("GenomeInfo_id", chromosomes), 
	FOREIGN KEY("GenomeInfo_id") REFERENCES "GenomeInfo" (id)
);
CREATE TABLE "LibrarySampleInfo_alternate_identifiers" (
	"LibrarySampleInfo_library_sample_name" TEXT, 
	alternate_identifiers TEXT, 
	PRIMARY KEY ("LibrarySampleInfo_library_sample_name", alternate_identifiers), 
	FOREIGN KEY("LibrarySampleInfo_library_sample_name") REFERENCES "LibrarySampleInfo" (library_sample_name)
);
CREATE TABLE "ProjectInfo_project_contributors" (
	"ProjectInfo_project_name" TEXT, 
	project_contributors TEXT, 
	PRIMARY KEY ("ProjectInfo_project_name", project_contributors), 
	FOREIGN KEY("ProjectInfo_project_name") REFERENCES "ProjectInfo" (project_name)
);
CREATE TABLE "SpecimenInfo_alternate_identifiers" (
	"SpecimenInfo_specimen_name" TEXT, 
	alternate_identifiers TEXT, 
	PRIMARY KEY ("SpecimenInfo_specimen_name", alternate_identifiers), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "SpecimenInfo_specimen_taxon_id" (
	"SpecimenInfo_specimen_name" TEXT, 
	specimen_taxon_id INTEGER NOT NULL, 
	PRIMARY KEY ("SpecimenInfo_specimen_name", specimen_taxon_id), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "SpecimenInfo_specimen_comments" (
	"SpecimenInfo_specimen_name" TEXT, 
	specimen_comments TEXT, 
	PRIMARY KEY ("SpecimenInfo_specimen_name", specimen_comments), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "SpecimenInfo_travel_out_six_month" (
	"SpecimenInfo_specimen_name" TEXT, 
	travel_out_six_month_id INTEGER, 
	PRIMARY KEY ("SpecimenInfo_specimen_name", travel_out_six_month_id), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name), 
	FOREIGN KEY(travel_out_six_month_id) REFERENCES "TravelInfo" (id)
);
CREATE TABLE "SpecimenInfo_drug_usage" (
	"SpecimenInfo_specimen_name" TEXT, 
	drug_usage TEXT, 
	PRIMARY KEY ("SpecimenInfo_specimen_name", drug_usage), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "SpecimenInfo_treatment_status" (
	"SpecimenInfo_specimen_name" TEXT, 
	treatment_status TEXT, 
	PRIMARY KEY ("SpecimenInfo_specimen_name", treatment_status), 
	FOREIGN KEY("SpecimenInfo_specimen_name") REFERENCES "SpecimenInfo" (specimen_name)
);
CREATE TABLE "DetectedMicrohaplotypesForTarget" (
	id INTEGER NOT NULL, 
	mhaps_target_id INTEGER NOT NULL, 
	"DetectedMicrohaplotypesForSample_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("DetectedMicrohaplotypesForSample_id") REFERENCES "DetectedMicrohaplotypesForSample" (id)
);
CREATE TABLE "ReadCountsByStageForTarget" (
	id INTEGER NOT NULL, 
	target_id INTEGER NOT NULL, 
	"ReadCountsByStageForLibrarySample_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("ReadCountsByStageForLibrarySample_id") REFERENCES "ReadCountsByStageForLibrarySample" (id)
);
CREATE TABLE "BioMethod_additional_argument" (
	"BioMethod_id" INTEGER, 
	additional_argument TEXT, 
	PRIMARY KEY ("BioMethod_id", additional_argument), 
	FOREIGN KEY("BioMethod_id") REFERENCES "BioMethod" (id)
);
CREATE TABLE "MicrohaplotypeForTarget" (
	id INTEGER NOT NULL, 
	mhap_id INTEGER NOT NULL, 
	reads INTEGER NOT NULL, 
	umis INTEGER, 
	"DetectedMicrohaplotypesForTarget_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("DetectedMicrohaplotypesForTarget_id") REFERENCES "DetectedMicrohaplotypesForTarget" (id)
);
CREATE TABLE "StageReadCounts" (
	id INTEGER NOT NULL, 
	reads INTEGER NOT NULL, 
	stage TEXT NOT NULL, 
	"ReadCountsByStageForTarget_id" INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY("ReadCountsByStageForTarget_id") REFERENCES "ReadCountsByStageForTarget" (id)
);