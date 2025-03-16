# Auto generated from portable_microhaplotype_object.yaml by pythongen.py version: 0.0.1
# Generation date: 2025-03-16T16:45:07
# Schema: portable-microhaplotype-object
#
# id: https://plasmogenepi.github.io/portable-microhaplotype-object
# description: A schema to define the minimum amount of data needed to export a microhaplotype calling pipeline analysis with associated metadata
# license: GNU GPL v3.0

import dataclasses
import re
from jsonasobj2 import JsonObj, as_dict
from typing import Optional, List, Union, Dict, ClassVar, Any
from dataclasses import dataclass
from datetime import date, datetime
from linkml_runtime.linkml_model.meta import EnumDefinition, PermissibleValue, PvFormulaOptions

from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.metamodelcore import empty_list, empty_dict, bnode
from linkml_runtime.utils.yamlutils import YAMLRoot, extended_str, extended_float, extended_int
from linkml_runtime.utils.dataclass_extensions_376 import dataclasses_init_fn_with_kwargs
from linkml_runtime.utils.formatutils import camelcase, underscore, sfx
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from rdflib import Namespace, URIRef
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.linkml_model.types import Double, Integer, String

metamodel_version = "1.7.0"
version = None

# Overwrite dataclasses _init_fn to add **kwargs in __init__
dataclasses._init_fn = dataclasses_init_fn_with_kwargs

# Namespaces
PATO = CurieNamespace('PATO', 'http://purl.obolibrary.org/obo/PATO_')
BIOLINK = CurieNamespace('biolink', 'https://w3id.org/biolink/')
EXAMPLE = CurieNamespace('example', 'https://example.org/')
LINKML = CurieNamespace('linkml', 'https://w3id.org/linkml/')
PORTABLE_MICROHAPLOTYPE_OBJECT = CurieNamespace('portable_microhaplotype_object', 'https://plasmogenepi.github.io/portable-microhaplotype-object/')
SCHEMA = CurieNamespace('schema', 'http://schema.org/')
DEFAULT_ = PORTABLE_MICROHAPLOTYPE_OBJECT


# Types

# Class references
class RepresentativeMicrohaplotypesForTargetTargetId(extended_int):
    pass


class ExperimentInfoExperimentSampleName(extended_str):
    pass


class SpecimenInfoSpecimenName(extended_str):
    pass


@dataclass
class MarkerOfInterest(YAMLRoot):
    """
    A specific genomic location of interest, e.g. drug resistance, or other phenotypical marker
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MarkerOfInterest"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MarkerOfInterest"
    class_name: ClassVar[str] = "MarkerOfInterest"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MarkerOfInterest

    marker_location: Optional[Union[dict, "GenomicLocation"]] = None
    associations: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.marker_location is not None and not isinstance(self.marker_location, GenomicLocation):
            self.marker_location = GenomicLocation(**as_dict(self.marker_location))

        if not isinstance(self.associations, list):
            self.associations = [self.associations] if self.associations is not None else []
        self.associations = [v if isinstance(v, str) else str(v) for v in self.associations]

        super().__post_init__(**kwargs)


@dataclass
class TargetInfo(YAMLRoot):
    """
    Information about a specific target within a genome
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["TargetInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:TargetInfo"
    class_name: ClassVar[str] = "TargetInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.TargetInfo

    target_name: str = None
    forward_primers: Union[Union[dict, "PrimerInfo"], List[Union[dict, "PrimerInfo"]]] = None
    reverse_primers: Union[Union[dict, "PrimerInfo"], List[Union[dict, "PrimerInfo"]]] = None
    gene_name: Optional[str] = None
    insert_location: Optional[Union[dict, "GenomicLocation"]] = None
    target_attributes: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_name):
            self.MissingRequiredField("target_name")
        if not isinstance(self.target_name, str):
            self.target_name = str(self.target_name)

        if self._is_empty(self.forward_primers):
            self.MissingRequiredField("forward_primers")
        self._normalize_inlined_as_dict(slot_name="forward_primers", slot_type=PrimerInfo, key_name="seq", keyed=False)

        if self._is_empty(self.reverse_primers):
            self.MissingRequiredField("reverse_primers")
        self._normalize_inlined_as_dict(slot_name="reverse_primers", slot_type=PrimerInfo, key_name="seq", keyed=False)

        if self.gene_name is not None and not isinstance(self.gene_name, str):
            self.gene_name = str(self.gene_name)

        if self.insert_location is not None and not isinstance(self.insert_location, GenomicLocation):
            self.insert_location = GenomicLocation(**as_dict(self.insert_location))

        if not isinstance(self.target_attributes, list):
            self.target_attributes = [self.target_attributes] if self.target_attributes is not None else []
        self.target_attributes = [v if isinstance(v, str) else str(v) for v in self.target_attributes]

        super().__post_init__(**kwargs)


@dataclass
class ReactionInfo(YAMLRoot):
    """
    information on a panel of targeted amplicon primer pairs
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReactionInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReactionInfo"
    class_name: ClassVar[str] = "ReactionInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReactionInfo

    panel_targets: Union[int, List[int]] = None
    reaction_name: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.panel_targets):
            self.MissingRequiredField("panel_targets")
        if not isinstance(self.panel_targets, list):
            self.panel_targets = [self.panel_targets] if self.panel_targets is not None else []
        self.panel_targets = [v if isinstance(v, int) else int(v) for v in self.panel_targets]

        if self._is_empty(self.reaction_name):
            self.MissingRequiredField("reaction_name")
        if not isinstance(self.reaction_name, str):
            self.reaction_name = str(self.reaction_name)

        super().__post_init__(**kwargs)


@dataclass
class PanelInfo(YAMLRoot):
    """
    information on a panel of targeted amplicon primer pairs
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PanelInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PanelInfo"
    class_name: ClassVar[str] = "PanelInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PanelInfo

    reactions: Union[Union[dict, ReactionInfo], List[Union[dict, ReactionInfo]]] = None
    panel_name: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.reactions):
            self.MissingRequiredField("reactions")
        self._normalize_inlined_as_dict(slot_name="reactions", slot_type=ReactionInfo, key_name="panel_targets", keyed=False)

        if self._is_empty(self.panel_name):
            self.MissingRequiredField("panel_name")
        if not isinstance(self.panel_name, str):
            self.panel_name = str(self.panel_name)

        super().__post_init__(**kwargs)


@dataclass
class RepresentativeMicrohaplotype(YAMLRoot):
    """
    the representative sequence for a microhaplotype, similar to a fast(a/q) format
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotype"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotype"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotype"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotype

    seq: str = None
    microhaplotype_name: Optional[str] = None
    quality: Optional[str] = None
    pseudocigar: Optional[str] = None
    masking: Optional[Union[Union[dict, "MaskingInfo"], List[Union[dict, "MaskingInfo"]]]] = empty_list()
    alt_annotations: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.seq):
            self.MissingRequiredField("seq")
        if not isinstance(self.seq, str):
            self.seq = str(self.seq)

        if self.microhaplotype_name is not None and not isinstance(self.microhaplotype_name, str):
            self.microhaplotype_name = str(self.microhaplotype_name)

        if self.quality is not None and not isinstance(self.quality, str):
            self.quality = str(self.quality)

        if self.pseudocigar is not None and not isinstance(self.pseudocigar, str):
            self.pseudocigar = str(self.pseudocigar)

        self._normalize_inlined_as_dict(slot_name="masking", slot_type=MaskingInfo, key_name="seq_start", keyed=False)

        if not isinstance(self.alt_annotations, list):
            self.alt_annotations = [self.alt_annotations] if self.alt_annotations is not None else []
        self.alt_annotations = [v if isinstance(v, str) else str(v) for v in self.alt_annotations]

        super().__post_init__(**kwargs)


@dataclass
class MaskingInfo(YAMLRoot):
    """
    information about a subsegment of the sequence that should be masked
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MaskingInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MaskingInfo"
    class_name: ClassVar[str] = "MaskingInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MaskingInfo

    seq_start: int = None
    seq_segment_size: int = None
    replacement_size: int = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.seq_start):
            self.MissingRequiredField("seq_start")
        if not isinstance(self.seq_start, int):
            self.seq_start = int(self.seq_start)

        if self._is_empty(self.seq_segment_size):
            self.MissingRequiredField("seq_segment_size")
        if not isinstance(self.seq_segment_size, int):
            self.seq_segment_size = int(self.seq_segment_size)

        if self._is_empty(self.replacement_size):
            self.MissingRequiredField("replacement_size")
        if not isinstance(self.replacement_size, int):
            self.replacement_size = int(self.replacement_size)

        super().__post_init__(**kwargs)


@dataclass
class RepresentativeMicrohaplotypes(YAMLRoot):
    """
    a collection of representative sequences for microhaplotypes for all targets
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotypes"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotypes"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotypes"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypes

    targets: Union[Dict[Union[int, RepresentativeMicrohaplotypesForTargetTargetId], Union[dict, "RepresentativeMicrohaplotypesForTarget"]], List[Union[dict, "RepresentativeMicrohaplotypesForTarget"]]] = empty_dict()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.targets):
            self.MissingRequiredField("targets")
        self._normalize_inlined_as_list(slot_name="targets", slot_type=RepresentativeMicrohaplotypesForTarget, key_name="target_id", keyed=True)

        super().__post_init__(**kwargs)


@dataclass
class RepresentativeMicrohaplotypesForTarget(YAMLRoot):
    """
    a list of the representative sequence for a microhaplotypes, similar to a fast(a/q) format
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotypesForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotypesForTarget"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotypesForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypesForTarget

    target_id: Union[int, RepresentativeMicrohaplotypesForTargetTargetId] = None
    microhaplotypes: Union[Union[dict, RepresentativeMicrohaplotype], List[Union[dict, RepresentativeMicrohaplotype]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, RepresentativeMicrohaplotypesForTargetTargetId):
            self.target_id = RepresentativeMicrohaplotypesForTargetTargetId(self.target_id)

        if self._is_empty(self.microhaplotypes):
            self.MissingRequiredField("microhaplotypes")
        self._normalize_inlined_as_dict(slot_name="microhaplotypes", slot_type=RepresentativeMicrohaplotype, key_name="seq", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class MicrohaplotypesDetected(YAMLRoot):
    """
    the microhaplotypes detected in a targeted amplicon analysis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MicrohaplotypesDetected"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MicrohaplotypesDetected"
    class_name: ClassVar[str] = "MicrohaplotypesDetected"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MicrohaplotypesDetected

    bioinformatics_run_id: int = None
    experiment_samples: Union[Union[dict, "MicrohaplotypesForSample"], List[Union[dict, "MicrohaplotypesForSample"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.bioinformatics_run_id):
            self.MissingRequiredField("bioinformatics_run_id")
        if not isinstance(self.bioinformatics_run_id, int):
            self.bioinformatics_run_id = int(self.bioinformatics_run_id)

        if self._is_empty(self.experiment_samples):
            self.MissingRequiredField("experiment_samples")
        if not isinstance(self.experiment_samples, list):
            self.experiment_samples = [self.experiment_samples] if self.experiment_samples is not None else []
        self.experiment_samples = [v if isinstance(v, MicrohaplotypesForSample) else MicrohaplotypesForSample(**as_dict(v)) for v in self.experiment_samples]

        super().__post_init__(**kwargs)


@dataclass
class GenomeInfo(YAMLRoot):
    """
    information on a genome
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["GenomeInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:GenomeInfo"
    class_name: ClassVar[str] = "GenomeInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.GenomeInfo

    name: str = None
    genome_version: str = None
    taxon_id: int = None
    url: str = None
    chromosomes: Optional[Union[str, List[str]]] = empty_list()
    gff_url: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.name):
            self.MissingRequiredField("name")
        if not isinstance(self.name, str):
            self.name = str(self.name)

        if self._is_empty(self.genome_version):
            self.MissingRequiredField("genome_version")
        if not isinstance(self.genome_version, str):
            self.genome_version = str(self.genome_version)

        if self._is_empty(self.taxon_id):
            self.MissingRequiredField("taxon_id")
        if not isinstance(self.taxon_id, int):
            self.taxon_id = int(self.taxon_id)

        if self._is_empty(self.url):
            self.MissingRequiredField("url")
        if not isinstance(self.url, str):
            self.url = str(self.url)

        if not isinstance(self.chromosomes, list):
            self.chromosomes = [self.chromosomes] if self.chromosomes is not None else []
        self.chromosomes = [v if isinstance(v, str) else str(v) for v in self.chromosomes]

        if self.gff_url is not None and not isinstance(self.gff_url, str):
            self.gff_url = str(self.gff_url)

        super().__post_init__(**kwargs)


@dataclass
class GenomicLocation(YAMLRoot):
    """
    information on the genomic location of specific sequence
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["GenomicLocation"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:GenomicLocation"
    class_name: ClassVar[str] = "GenomicLocation"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.GenomicLocation

    genome_id: int = None
    chrom: str = None
    start: int = None
    end: int = None
    strand: Optional[str] = None
    ref_seq: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.genome_id):
            self.MissingRequiredField("genome_id")
        if not isinstance(self.genome_id, int):
            self.genome_id = int(self.genome_id)

        if self._is_empty(self.chrom):
            self.MissingRequiredField("chrom")
        if not isinstance(self.chrom, str):
            self.chrom = str(self.chrom)

        if self._is_empty(self.start):
            self.MissingRequiredField("start")
        if not isinstance(self.start, int):
            self.start = int(self.start)

        if self._is_empty(self.end):
            self.MissingRequiredField("end")
        if not isinstance(self.end, int):
            self.end = int(self.end)

        if self.strand is not None and not isinstance(self.strand, str):
            self.strand = str(self.strand)

        if self.ref_seq is not None and not isinstance(self.ref_seq, str):
            self.ref_seq = str(self.ref_seq)

        super().__post_init__(**kwargs)


@dataclass
class PrimerInfo(YAMLRoot):
    """
    information on a primer sequence
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PrimerInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PrimerInfo"
    class_name: ClassVar[str] = "PrimerInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PrimerInfo

    seq: str = None
    location: Optional[Union[dict, GenomicLocation]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.seq):
            self.MissingRequiredField("seq")
        if not isinstance(self.seq, str):
            self.seq = str(self.seq)

        if self.location is not None and not isinstance(self.location, GenomicLocation):
            self.location = GenomicLocation(**as_dict(self.location))

        super().__post_init__(**kwargs)


@dataclass
class MicrohaplotypesForSample(YAMLRoot):
    """
    Microhaplotypes detected for a sample for all targets
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MicrohaplotypesForSample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MicrohaplotypesForSample"
    class_name: ClassVar[str] = "MicrohaplotypesForSample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MicrohaplotypesForSample

    experiment_sample_id: int = None
    target_results: Union[Union[dict, "MicrohaplotypesForTarget"], List[Union[dict, "MicrohaplotypesForTarget"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_sample_id):
            self.MissingRequiredField("experiment_sample_id")
        if not isinstance(self.experiment_sample_id, int):
            self.experiment_sample_id = int(self.experiment_sample_id)

        if self._is_empty(self.target_results):
            self.MissingRequiredField("target_results")
        if not isinstance(self.target_results, list):
            self.target_results = [self.target_results] if self.target_results is not None else []
        self.target_results = [v if isinstance(v, MicrohaplotypesForTarget) else MicrohaplotypesForTarget(**as_dict(v)) for v in self.target_results]

        super().__post_init__(**kwargs)


@dataclass
class MicrohaplotypeForTarget(YAMLRoot):
    """
    Microhaplotype detected for a specific target
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MicrohaplotypeForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MicrohaplotypeForTarget"
    class_name: ClassVar[str] = "MicrohaplotypeForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MicrohaplotypeForTarget

    mhap_id: int = None
    reads: int = None
    umis: Optional[int] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.mhap_id):
            self.MissingRequiredField("mhap_id")
        if not isinstance(self.mhap_id, int):
            self.mhap_id = int(self.mhap_id)

        if self._is_empty(self.reads):
            self.MissingRequiredField("reads")
        if not isinstance(self.reads, int):
            self.reads = int(self.reads)

        if self.umis is not None and not isinstance(self.umis, int):
            self.umis = int(self.umis)

        super().__post_init__(**kwargs)


@dataclass
class MicrohaplotypesForTarget(YAMLRoot):
    """
    Microhaplotypes detected for a specific target
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MicrohaplotypesForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MicrohaplotypesForTarget"
    class_name: ClassVar[str] = "MicrohaplotypesForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MicrohaplotypesForTarget

    mhaps_target_id: int = None
    haps: Union[Union[dict, MicrohaplotypeForTarget], List[Union[dict, MicrohaplotypeForTarget]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.mhaps_target_id):
            self.MissingRequiredField("mhaps_target_id")
        if not isinstance(self.mhaps_target_id, int):
            self.mhaps_target_id = int(self.mhaps_target_id)

        if self._is_empty(self.haps):
            self.MissingRequiredField("haps")
        if not isinstance(self.haps, list):
            self.haps = [self.haps] if self.haps is not None else []
        self.haps = [v if isinstance(v, MicrohaplotypeForTarget) else MicrohaplotypeForTarget(**as_dict(v)) for v in self.haps]

        super().__post_init__(**kwargs)


@dataclass
class BioinformaticsMethodInfo(YAMLRoot):
    """
    the targeted amplicon bioinformatics pipeline
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["BioinformaticsMethodInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:BioinformaticsMethodInfo"
    class_name: ClassVar[str] = "BioinformaticsMethodInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.BioinformaticsMethodInfo

    demultiplexing_method: Union[dict, "BioMethod"] = None
    denoising_method: Union[dict, "BioMethod"] = None
    additional_methods: Optional[Union[Union[dict, "BioMethod"], List[Union[dict, "BioMethod"]]]] = empty_list()
    bioinformatics_method_name: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.demultiplexing_method):
            self.MissingRequiredField("demultiplexing_method")
        if not isinstance(self.demultiplexing_method, BioMethod):
            self.demultiplexing_method = BioMethod(**as_dict(self.demultiplexing_method))

        if self._is_empty(self.denoising_method):
            self.MissingRequiredField("denoising_method")
        if not isinstance(self.denoising_method, BioMethod):
            self.denoising_method = BioMethod(**as_dict(self.denoising_method))

        if not isinstance(self.additional_methods, list):
            self.additional_methods = [self.additional_methods] if self.additional_methods is not None else []
        self.additional_methods = [v if isinstance(v, BioMethod) else BioMethod(**as_dict(v)) for v in self.additional_methods]

        if self.bioinformatics_method_name is not None and not isinstance(self.bioinformatics_method_name, str):
            self.bioinformatics_method_name = str(self.bioinformatics_method_name)

        super().__post_init__(**kwargs)


@dataclass
class BioMethod(YAMLRoot):
    """
    methodology description of a portion of a bioinformatics pipeline
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["BioMethod"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:BioMethod"
    class_name: ClassVar[str] = "BioMethod"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.BioMethod

    program_version: str = None
    program: str = None
    program_description: Optional[str] = None
    additional_argument: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.program_version):
            self.MissingRequiredField("program_version")
        if not isinstance(self.program_version, str):
            self.program_version = str(self.program_version)

        if self._is_empty(self.program):
            self.MissingRequiredField("program")
        if not isinstance(self.program, str):
            self.program = str(self.program)

        if self.program_description is not None and not isinstance(self.program_description, str):
            self.program_description = str(self.program_description)

        if not isinstance(self.additional_argument, list):
            self.additional_argument = [self.additional_argument] if self.additional_argument is not None else []
        self.additional_argument = [v if isinstance(v, str) else str(v) for v in self.additional_argument]

        super().__post_init__(**kwargs)


@dataclass
class ExperimentInfo(YAMLRoot):
    """
    Information about a specific amplification and sequencing of a specimen
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ExperimentInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ExperimentInfo"
    class_name: ClassVar[str] = "ExperimentInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo

    experiment_sample_name: Union[str, ExperimentInfoExperimentSampleName] = None
    sequencing_info_id: int = None
    specimen_id: int = None
    panel_id: int = None
    plate_name: Optional[str] = None
    plate_row: Optional[str] = None
    plate_col: Optional[int] = None
    accession: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_sample_name):
            self.MissingRequiredField("experiment_sample_name")
        if not isinstance(self.experiment_sample_name, ExperimentInfoExperimentSampleName):
            self.experiment_sample_name = ExperimentInfoExperimentSampleName(self.experiment_sample_name)

        if self._is_empty(self.sequencing_info_id):
            self.MissingRequiredField("sequencing_info_id")
        if not isinstance(self.sequencing_info_id, int):
            self.sequencing_info_id = int(self.sequencing_info_id)

        if self._is_empty(self.specimen_id):
            self.MissingRequiredField("specimen_id")
        if not isinstance(self.specimen_id, int):
            self.specimen_id = int(self.specimen_id)

        if self._is_empty(self.panel_id):
            self.MissingRequiredField("panel_id")
        if not isinstance(self.panel_id, int):
            self.panel_id = int(self.panel_id)

        if self.plate_name is not None and not isinstance(self.plate_name, str):
            self.plate_name = str(self.plate_name)

        if self.plate_row is not None and not isinstance(self.plate_row, str):
            self.plate_row = str(self.plate_row)

        if self.plate_col is not None and not isinstance(self.plate_col, int):
            self.plate_col = int(self.plate_col)

        if self.accession is not None and not isinstance(self.accession, str):
            self.accession = str(self.accession)

        super().__post_init__(**kwargs)


@dataclass
class SequencingInfo(YAMLRoot):
    """
    Information on sequencing info
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["SequencingInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:SequencingInfo"
    class_name: ClassVar[str] = "SequencingInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.SequencingInfo

    sequencing_info_name: str = None
    seq_instrument: str = None
    seq_date: str = None
    nucl_acid_ext: Optional[str] = None
    nucl_acid_amp: Optional[str] = None
    nucl_acid_ext_date: Optional[str] = None
    nucl_acid_amp_date: Optional[str] = None
    pcr_cond: Optional[str] = None
    lib_screen: Optional[str] = None
    lib_layout: Optional[str] = None
    lib_kit: Optional[str] = None
    seq_center: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.sequencing_info_name):
            self.MissingRequiredField("sequencing_info_name")
        if not isinstance(self.sequencing_info_name, str):
            self.sequencing_info_name = str(self.sequencing_info_name)

        if self._is_empty(self.seq_instrument):
            self.MissingRequiredField("seq_instrument")
        if not isinstance(self.seq_instrument, str):
            self.seq_instrument = str(self.seq_instrument)

        if self._is_empty(self.seq_date):
            self.MissingRequiredField("seq_date")
        if not isinstance(self.seq_date, str):
            self.seq_date = str(self.seq_date)

        if self.nucl_acid_ext is not None and not isinstance(self.nucl_acid_ext, str):
            self.nucl_acid_ext = str(self.nucl_acid_ext)

        if self.nucl_acid_amp is not None and not isinstance(self.nucl_acid_amp, str):
            self.nucl_acid_amp = str(self.nucl_acid_amp)

        if self.nucl_acid_ext_date is not None and not isinstance(self.nucl_acid_ext_date, str):
            self.nucl_acid_ext_date = str(self.nucl_acid_ext_date)

        if self.nucl_acid_amp_date is not None and not isinstance(self.nucl_acid_amp_date, str):
            self.nucl_acid_amp_date = str(self.nucl_acid_amp_date)

        if self.pcr_cond is not None and not isinstance(self.pcr_cond, str):
            self.pcr_cond = str(self.pcr_cond)

        if self.lib_screen is not None and not isinstance(self.lib_screen, str):
            self.lib_screen = str(self.lib_screen)

        if self.lib_layout is not None and not isinstance(self.lib_layout, str):
            self.lib_layout = str(self.lib_layout)

        if self.lib_kit is not None and not isinstance(self.lib_kit, str):
            self.lib_kit = str(self.lib_kit)

        if self.seq_center is not None and not isinstance(self.seq_center, str):
            self.seq_center = str(self.seq_center)

        super().__post_init__(**kwargs)


@dataclass
class ParasiteDensity(YAMLRoot):
    """
    method and value of determined parasite density
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ParasiteDensity"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ParasiteDensity"
    class_name: ClassVar[str] = "ParasiteDensity"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ParasiteDensity

    method: str = None
    density: float = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.method):
            self.MissingRequiredField("method")
        if not isinstance(self.method, str):
            self.method = str(self.method)

        if self._is_empty(self.density):
            self.MissingRequiredField("density")
        if not isinstance(self.density, float):
            self.density = float(self.density)

        super().__post_init__(**kwargs)


@dataclass
class SpecimenInfo(YAMLRoot):
    """
    Information on specimen info
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["SpecimenInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:SpecimenInfo"
    class_name: ClassVar[str] = "SpecimenInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.SpecimenInfo

    specimen_name: Union[str, SpecimenInfoSpecimenName] = None
    samp_taxon_id: int = None
    collection_date: str = None
    collection_country: str = None
    collector: str = None
    samp_store_loc: str = None
    samp_collect_device: str = None
    project_name: str = None
    plate_name: Optional[str] = None
    plate_row: Optional[str] = None
    plate_col: Optional[int] = None
    individual_id: Optional[int] = None
    host_taxon_id: Optional[int] = None
    alternate_identifiers: Optional[Union[str, List[str]]] = empty_list()
    specimen_sex: Optional[str] = None
    parasite_density_info: Optional[Union[Union[dict, ParasiteDensity], List[Union[dict, ParasiteDensity]]]] = empty_list()
    date_of_birth: Optional[str] = None
    geo_admin1: Optional[str] = None
    geo_admin2: Optional[str] = None
    geo_admin3: Optional[str] = None
    lat_lon: Optional[str] = None
    sample_comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.specimen_name):
            self.MissingRequiredField("specimen_name")
        if not isinstance(self.specimen_name, SpecimenInfoSpecimenName):
            self.specimen_name = SpecimenInfoSpecimenName(self.specimen_name)

        if self._is_empty(self.samp_taxon_id):
            self.MissingRequiredField("samp_taxon_id")
        if not isinstance(self.samp_taxon_id, int):
            self.samp_taxon_id = int(self.samp_taxon_id)

        if self._is_empty(self.collection_date):
            self.MissingRequiredField("collection_date")
        if not isinstance(self.collection_date, str):
            self.collection_date = str(self.collection_date)

        if self._is_empty(self.collection_country):
            self.MissingRequiredField("collection_country")
        if not isinstance(self.collection_country, str):
            self.collection_country = str(self.collection_country)

        if self._is_empty(self.collector):
            self.MissingRequiredField("collector")
        if not isinstance(self.collector, str):
            self.collector = str(self.collector)

        if self._is_empty(self.samp_store_loc):
            self.MissingRequiredField("samp_store_loc")
        if not isinstance(self.samp_store_loc, str):
            self.samp_store_loc = str(self.samp_store_loc)

        if self._is_empty(self.samp_collect_device):
            self.MissingRequiredField("samp_collect_device")
        if not isinstance(self.samp_collect_device, str):
            self.samp_collect_device = str(self.samp_collect_device)

        if self._is_empty(self.project_name):
            self.MissingRequiredField("project_name")
        if not isinstance(self.project_name, str):
            self.project_name = str(self.project_name)

        if self.plate_name is not None and not isinstance(self.plate_name, str):
            self.plate_name = str(self.plate_name)

        if self.plate_row is not None and not isinstance(self.plate_row, str):
            self.plate_row = str(self.plate_row)

        if self.plate_col is not None and not isinstance(self.plate_col, int):
            self.plate_col = int(self.plate_col)

        if self.individual_id is not None and not isinstance(self.individual_id, int):
            self.individual_id = int(self.individual_id)

        if self.host_taxon_id is not None and not isinstance(self.host_taxon_id, int):
            self.host_taxon_id = int(self.host_taxon_id)

        if not isinstance(self.alternate_identifiers, list):
            self.alternate_identifiers = [self.alternate_identifiers] if self.alternate_identifiers is not None else []
        self.alternate_identifiers = [v if isinstance(v, str) else str(v) for v in self.alternate_identifiers]

        if self.specimen_sex is not None and not isinstance(self.specimen_sex, str):
            self.specimen_sex = str(self.specimen_sex)

        if not isinstance(self.parasite_density_info, list):
            self.parasite_density_info = [self.parasite_density_info] if self.parasite_density_info is not None else []
        self.parasite_density_info = [v if isinstance(v, ParasiteDensity) else ParasiteDensity(**as_dict(v)) for v in self.parasite_density_info]

        if self.date_of_birth is not None and not isinstance(self.date_of_birth, str):
            self.date_of_birth = str(self.date_of_birth)

        if self.geo_admin1 is not None and not isinstance(self.geo_admin1, str):
            self.geo_admin1 = str(self.geo_admin1)

        if self.geo_admin2 is not None and not isinstance(self.geo_admin2, str):
            self.geo_admin2 = str(self.geo_admin2)

        if self.geo_admin3 is not None and not isinstance(self.geo_admin3, str):
            self.geo_admin3 = str(self.geo_admin3)

        if self.lat_lon is not None and not isinstance(self.lat_lon, str):
            self.lat_lon = str(self.lat_lon)

        if self.sample_comments is not None and not isinstance(self.sample_comments, str):
            self.sample_comments = str(self.sample_comments)

        super().__post_init__(**kwargs)


@dataclass
class BioinformaticsRunInfo(YAMLRoot):
    """
    Information about the pipeline run that generated some of the microhaplotype detected and reads_by_stage
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["BioinformaticsRunInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:BioinformaticsRunInfo"
    class_name: ClassVar[str] = "BioinformaticsRunInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.BioinformaticsRunInfo

    bioinformatics_methods_id: int = None
    run_date: str = None
    bioinformatics_run_name: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.bioinformatics_methods_id):
            self.MissingRequiredField("bioinformatics_methods_id")
        if not isinstance(self.bioinformatics_methods_id, int):
            self.bioinformatics_methods_id = int(self.bioinformatics_methods_id)

        if self._is_empty(self.run_date):
            self.MissingRequiredField("run_date")
        if not isinstance(self.run_date, str):
            self.run_date = str(self.run_date)

        if self.bioinformatics_run_name is not None and not isinstance(self.bioinformatics_run_name, str):
            self.bioinformatics_run_name = str(self.bioinformatics_run_name)

        super().__post_init__(**kwargs)


@dataclass
class PmoGenerationMethod(YAMLRoot):
    """
    Information about how a PMO was generated
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PmoGenerationMethod"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PmoGenerationMethod"
    class_name: ClassVar[str] = "PmoGenerationMethod"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PmoGenerationMethod

    program_version: str = None
    program_name: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.program_version):
            self.MissingRequiredField("program_version")
        if not isinstance(self.program_version, str):
            self.program_version = str(self.program_version)

        if self._is_empty(self.program_name):
            self.MissingRequiredField("program_name")
        if not isinstance(self.program_name, str):
            self.program_name = str(self.program_name)

        super().__post_init__(**kwargs)


@dataclass
class PmoHeader(YAMLRoot):
    """
    Information on the PMO file
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PmoHeader"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PmoHeader"
    class_name: ClassVar[str] = "PmoHeader"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PmoHeader

    pmo_version: str = None
    creation_date: Optional[str] = None
    generation_method: Optional[Union[dict, PmoGenerationMethod]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.pmo_version):
            self.MissingRequiredField("pmo_version")
        if not isinstance(self.pmo_version, str):
            self.pmo_version = str(self.pmo_version)

        if self.creation_date is not None and not isinstance(self.creation_date, str):
            self.creation_date = str(self.creation_date)

        if self.generation_method is not None and not isinstance(self.generation_method, PmoGenerationMethod):
            self.generation_method = PmoGenerationMethod(**as_dict(self.generation_method))

        super().__post_init__(**kwargs)


@dataclass
class StageReadCounts(YAMLRoot):
    """
    Information on the reads counts at several stages
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["StageReadCounts"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:StageReadCounts"
    class_name: ClassVar[str] = "StageReadCounts"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.StageReadCounts

    read_count: int = None
    stage: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.read_count):
            self.MissingRequiredField("read_count")
        if not isinstance(self.read_count, int):
            self.read_count = int(self.read_count)

        if self._is_empty(self.stage):
            self.MissingRequiredField("stage")
        if not isinstance(self.stage, str):
            self.stage = str(self.stage)

        super().__post_init__(**kwargs)


@dataclass
class ReadCountsByStageForTarget(YAMLRoot):
    """
    Information on the reads counts at several stages of a pipeline for a target
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReadCountsByStageForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReadCountsByStageForTarget"
    class_name: ClassVar[str] = "ReadCountsByStageForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReadCountsByStageForTarget

    target_id: int = None
    stages: Union[Union[dict, StageReadCounts], List[Union[dict, StageReadCounts]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, int):
            self.target_id = int(self.target_id)

        if self._is_empty(self.stages):
            self.MissingRequiredField("stages")
        if not isinstance(self.stages, list):
            self.stages = [self.stages] if self.stages is not None else []
        self.stages = [v if isinstance(v, StageReadCounts) else StageReadCounts(**as_dict(v)) for v in self.stages]

        super().__post_init__(**kwargs)


@dataclass
class ReadCountsByStageForExperimentalSample(YAMLRoot):
    """
    Information on the reads counts at several stages of a pipeline for a experimental_sample
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReadCountsByStageForExperimentalSample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReadCountsByStageForExperimentalSample"
    class_name: ClassVar[str] = "ReadCountsByStageForExperimentalSample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReadCountsByStageForExperimentalSample

    experiment_sample_id: int = None
    total_raw_count: int = None
    read_counts_for_targets: Optional[Union[Union[dict, ReadCountsByStageForTarget], List[Union[dict, ReadCountsByStageForTarget]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_sample_id):
            self.MissingRequiredField("experiment_sample_id")
        if not isinstance(self.experiment_sample_id, int):
            self.experiment_sample_id = int(self.experiment_sample_id)

        if self._is_empty(self.total_raw_count):
            self.MissingRequiredField("total_raw_count")
        if not isinstance(self.total_raw_count, int):
            self.total_raw_count = int(self.total_raw_count)

        if not isinstance(self.read_counts_for_targets, list):
            self.read_counts_for_targets = [self.read_counts_for_targets] if self.read_counts_for_targets is not None else []
        self.read_counts_for_targets = [v if isinstance(v, ReadCountsByStageForTarget) else ReadCountsByStageForTarget(**as_dict(v)) for v in self.read_counts_for_targets]

        super().__post_init__(**kwargs)


@dataclass
class ReadCountsByStage(YAMLRoot):
    """
    Information on the reads counts at several stages of a pipeline
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReadCountsByStage"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReadCountsByStage"
    class_name: ClassVar[str] = "ReadCountsByStage"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReadCountsByStage

    bioinformatics_run_id: int = None
    read_counts_by_experimental_sample_by_stage: Union[Union[dict, ReadCountsByStageForExperimentalSample], List[Union[dict, ReadCountsByStageForExperimentalSample]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.bioinformatics_run_id):
            self.MissingRequiredField("bioinformatics_run_id")
        if not isinstance(self.bioinformatics_run_id, int):
            self.bioinformatics_run_id = int(self.bioinformatics_run_id)

        if self._is_empty(self.read_counts_by_experimental_sample_by_stage):
            self.MissingRequiredField("read_counts_by_experimental_sample_by_stage")
        if not isinstance(self.read_counts_by_experimental_sample_by_stage, list):
            self.read_counts_by_experimental_sample_by_stage = [self.read_counts_by_experimental_sample_by_stage] if self.read_counts_by_experimental_sample_by_stage is not None else []
        self.read_counts_by_experimental_sample_by_stage = [v if isinstance(v, ReadCountsByStageForExperimentalSample) else ReadCountsByStageForExperimentalSample(**as_dict(v)) for v in self.read_counts_by_experimental_sample_by_stage]

        super().__post_init__(**kwargs)


@dataclass
class PortableMicrohaplotypeObject(YAMLRoot):
    """
    Information on final results from a targeted amplicon analysis
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PortableMicrohaplotypeObject"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PortableMicrohaplotypeObject"
    class_name: ClassVar[str] = "PortableMicrohaplotypeObject"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PortableMicrohaplotypeObject

    experiment_info: Union[Dict[Union[str, ExperimentInfoExperimentSampleName], Union[dict, ExperimentInfo]], List[Union[dict, ExperimentInfo]]] = empty_dict()
    specimen_info: Union[Dict[Union[str, SpecimenInfoSpecimenName], Union[dict, SpecimenInfo]], List[Union[dict, SpecimenInfo]]] = empty_dict()
    sequencing_info: Union[Union[dict, SequencingInfo], List[Union[dict, SequencingInfo]]] = None
    panel_info: Union[Union[dict, PanelInfo], List[Union[dict, PanelInfo]]] = None
    target_info: Union[Union[dict, TargetInfo], List[Union[dict, TargetInfo]]] = None
    targeted_genomes: Union[Union[dict, GenomeInfo], List[Union[dict, GenomeInfo]]] = None
    microhaplotypes_info: Union[dict, RepresentativeMicrohaplotypes] = None
    bioinformatics_methods_info: Union[Union[dict, BioinformaticsMethodInfo], List[Union[dict, BioinformaticsMethodInfo]]] = None
    bioinformatics_run_info: Union[Union[dict, BioinformaticsRunInfo], List[Union[dict, BioinformaticsRunInfo]]] = None
    microhaplotypes_detected: Union[Union[dict, MicrohaplotypesDetected], List[Union[dict, MicrohaplotypesDetected]]] = None
    pmo_header: Union[dict, PmoHeader] = None
    read_counts_by_stage: Optional[Union[Union[dict, ReadCountsByStage], List[Union[dict, ReadCountsByStage]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_info):
            self.MissingRequiredField("experiment_info")
        self._normalize_inlined_as_list(slot_name="experiment_info", slot_type=ExperimentInfo, key_name="experiment_sample_name", keyed=True)

        if self._is_empty(self.specimen_info):
            self.MissingRequiredField("specimen_info")
        self._normalize_inlined_as_list(slot_name="specimen_info", slot_type=SpecimenInfo, key_name="specimen_name", keyed=True)

        if self._is_empty(self.sequencing_info):
            self.MissingRequiredField("sequencing_info")
        if not isinstance(self.sequencing_info, list):
            self.sequencing_info = [self.sequencing_info] if self.sequencing_info is not None else []
        self.sequencing_info = [v if isinstance(v, SequencingInfo) else SequencingInfo(**as_dict(v)) for v in self.sequencing_info]

        if self._is_empty(self.panel_info):
            self.MissingRequiredField("panel_info")
        if not isinstance(self.panel_info, list):
            self.panel_info = [self.panel_info] if self.panel_info is not None else []
        self.panel_info = [v if isinstance(v, PanelInfo) else PanelInfo(**as_dict(v)) for v in self.panel_info]

        if self._is_empty(self.target_info):
            self.MissingRequiredField("target_info")
        if not isinstance(self.target_info, list):
            self.target_info = [self.target_info] if self.target_info is not None else []
        self.target_info = [v if isinstance(v, TargetInfo) else TargetInfo(**as_dict(v)) for v in self.target_info]

        if self._is_empty(self.targeted_genomes):
            self.MissingRequiredField("targeted_genomes")
        if not isinstance(self.targeted_genomes, list):
            self.targeted_genomes = [self.targeted_genomes] if self.targeted_genomes is not None else []
        self.targeted_genomes = [v if isinstance(v, GenomeInfo) else GenomeInfo(**as_dict(v)) for v in self.targeted_genomes]

        if self._is_empty(self.microhaplotypes_info):
            self.MissingRequiredField("microhaplotypes_info")
        if not isinstance(self.microhaplotypes_info, RepresentativeMicrohaplotypes):
            self.microhaplotypes_info = RepresentativeMicrohaplotypes(**as_dict(self.microhaplotypes_info))

        if self._is_empty(self.bioinformatics_methods_info):
            self.MissingRequiredField("bioinformatics_methods_info")
        if not isinstance(self.bioinformatics_methods_info, list):
            self.bioinformatics_methods_info = [self.bioinformatics_methods_info] if self.bioinformatics_methods_info is not None else []
        self.bioinformatics_methods_info = [v if isinstance(v, BioinformaticsMethodInfo) else BioinformaticsMethodInfo(**as_dict(v)) for v in self.bioinformatics_methods_info]

        if self._is_empty(self.bioinformatics_run_info):
            self.MissingRequiredField("bioinformatics_run_info")
        if not isinstance(self.bioinformatics_run_info, list):
            self.bioinformatics_run_info = [self.bioinformatics_run_info] if self.bioinformatics_run_info is not None else []
        self.bioinformatics_run_info = [v if isinstance(v, BioinformaticsRunInfo) else BioinformaticsRunInfo(**as_dict(v)) for v in self.bioinformatics_run_info]

        if self._is_empty(self.microhaplotypes_detected):
            self.MissingRequiredField("microhaplotypes_detected")
        if not isinstance(self.microhaplotypes_detected, list):
            self.microhaplotypes_detected = [self.microhaplotypes_detected] if self.microhaplotypes_detected is not None else []
        self.microhaplotypes_detected = [v if isinstance(v, MicrohaplotypesDetected) else MicrohaplotypesDetected(**as_dict(v)) for v in self.microhaplotypes_detected]

        if self._is_empty(self.pmo_header):
            self.MissingRequiredField("pmo_header")
        if not isinstance(self.pmo_header, PmoHeader):
            self.pmo_header = PmoHeader(**as_dict(self.pmo_header))

        if not isinstance(self.read_counts_by_stage, list):
            self.read_counts_by_stage = [self.read_counts_by_stage] if self.read_counts_by_stage is not None else []
        self.read_counts_by_stage = [v if isinstance(v, ReadCountsByStage) else ReadCountsByStage(**as_dict(v)) for v in self.read_counts_by_stage]

        super().__post_init__(**kwargs)


# Enumerations


# Slots
class slots:
    pass

slots.experiment_sample_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_sample_id, name="experiment_sample_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_sample_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_sample_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.sequencing_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, name="sequencing_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.target_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_id, name="target_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.panel_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_id, name="panel_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.mhaps_target_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.mhaps_target_id, name="mhaps_target_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('mhaps_target_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.mhaps_target_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.mhap_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.mhap_id, name="mhap_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('mhap_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.mhap_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.bioinformatics_methods_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_methods_id, name="bioinformatics_methods_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_methods_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_methods_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.bioinformatics_run_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_run_id, name="bioinformatics_run_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_run_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_run_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.genome_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genome_id, name="genome_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('genome_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genome_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.plate_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, name="plate_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.plate_row = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, name="plate_row", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_row'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z]$'))

slots.plate_col = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, name="plate_col", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_col'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]+$'))

slots.seq = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq, name="seq", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq, domain=None, range=str,
                   pattern=re.compile(r'^[A-z]+$'))

slots.program_version = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_version, name="program_version", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program_version'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_version, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.markerOfInterest__marker_location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.marker_location, name="markerOfInterest__marker_location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('marker_location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.markerOfInterest__marker_location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.markerOfInterest__associations = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.associations, name="markerOfInterest__associations", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('associations'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.markerOfInterest__associations, domain=None, range=Optional[Union[str, List[str]]],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.targetInfo__target_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_name, name="targetInfo__target_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__target_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.targetInfo__gene_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gene_name, name="targetInfo__gene_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gene_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__gene_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.targetInfo__insert_location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.insert_location, name="targetInfo__insert_location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('insert_location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__insert_location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.targetInfo__forward_primers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.forward_primers, name="targetInfo__forward_primers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('forward_primers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__forward_primers, domain=None, range=Union[Union[dict, PrimerInfo], List[Union[dict, PrimerInfo]]])

slots.targetInfo__reverse_primers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reverse_primers, name="targetInfo__reverse_primers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reverse_primers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__reverse_primers, domain=None, range=Union[Union[dict, PrimerInfo], List[Union[dict, PrimerInfo]]])

slots.targetInfo__target_attributes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_attributes, name="targetInfo__target_attributes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_attributes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__target_attributes, domain=None, range=Optional[Union[str, List[str]]])

slots.reactionInfo__panel_targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_targets, name="reactionInfo__panel_targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reactionInfo__panel_targets, domain=None, range=Union[int, List[int]])

slots.reactionInfo__reaction_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reaction_name, name="reactionInfo__reaction_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reaction_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reactionInfo__reaction_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.panelInfo__reactions = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reactions, name="panelInfo__reactions", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reactions'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__reactions, domain=None, range=Union[Union[dict, ReactionInfo], List[Union[dict, ReactionInfo]]])

slots.panelInfo__panel_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_name, name="panelInfo__panel_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__panel_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.representativeMicrohaplotype__microhaplotype_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotype_name, name="representativeMicrohaplotype__microhaplotype_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotype_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__microhaplotype_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.representativeMicrohaplotype__quality = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.quality, name="representativeMicrohaplotype__quality", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('quality'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__quality, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.representativeMicrohaplotype__pseudocigar = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar, name="representativeMicrohaplotype__pseudocigar", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pseudocigar'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__pseudocigar, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.representativeMicrohaplotype__masking = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.masking, name="representativeMicrohaplotype__masking", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('masking'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__masking, domain=None, range=Optional[Union[Union[dict, MaskingInfo], List[Union[dict, MaskingInfo]]]])

slots.representativeMicrohaplotype__alt_annotations = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alt_annotations, name="representativeMicrohaplotype__alt_annotations", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alt_annotations'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__alt_annotations, domain=None, range=Optional[Union[str, List[str]]])

slots.maskingInfo__seq_start = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_start, name="maskingInfo__seq_start", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_start'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__seq_start, domain=None, range=int)

slots.maskingInfo__seq_segment_size = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_segment_size, name="maskingInfo__seq_segment_size", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_segment_size'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__seq_segment_size, domain=None, range=int)

slots.maskingInfo__replacement_size = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.replacement_size, name="maskingInfo__replacement_size", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('replacement_size'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__replacement_size, domain=None, range=int)

slots.representativeMicrohaplotypes__targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targets, name="representativeMicrohaplotypes__targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypes__targets, domain=None, range=Union[Dict[Union[int, RepresentativeMicrohaplotypesForTargetTargetId], Union[dict, RepresentativeMicrohaplotypesForTarget]], List[Union[dict, RepresentativeMicrohaplotypesForTarget]]])

slots.representativeMicrohaplotypesForTarget__microhaplotypes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes, name="representativeMicrohaplotypesForTarget__microhaplotypes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypesForTarget__microhaplotypes, domain=None, range=Union[Union[dict, RepresentativeMicrohaplotype], List[Union[dict, RepresentativeMicrohaplotype]]])

slots.microhaplotypesDetected__experiment_samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_samples, name="microhaplotypesDetected__experiment_samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesDetected__experiment_samples, domain=None, range=Union[Union[dict, MicrohaplotypesForSample], List[Union[dict, MicrohaplotypesForSample]]])

slots.genomeInfo__name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.name, name="genomeInfo__name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.genomeInfo__genome_version = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genome_version, name="genomeInfo__genome_version", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('genome_version'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__genome_version, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.genomeInfo__taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.taxon_id, name="genomeInfo__taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__taxon_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]$'))

slots.genomeInfo__url = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.url, name="genomeInfo__url", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('url'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__url, domain=None, range=str,
                   pattern=re.compile(r'^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$'))

slots.genomeInfo__chromosomes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.chromosomes, name="genomeInfo__chromosomes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('chromosomes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__chromosomes, domain=None, range=Optional[Union[str, List[str]]],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.genomeInfo__gff_url = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gff_url, name="genomeInfo__gff_url", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gff_url'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__gff_url, domain=None, range=Optional[str],
                   pattern=re.compile(r'^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$'))

slots.genomicLocation__chrom = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.chrom, name="genomicLocation__chrom", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('chrom'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__chrom, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.genomicLocation__start = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.start, name="genomicLocation__start", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('start'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__start, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.genomicLocation__end = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.end, name="genomicLocation__end", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('end'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__end, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.genomicLocation__strand = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.strand, name="genomicLocation__strand", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('strand'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__strand, domain=None, range=Optional[str],
                   pattern=re.compile(r'[+-]'))

slots.genomicLocation__ref_seq = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ref_seq, name="genomicLocation__ref_seq", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('ref_seq'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__ref_seq, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-]+$'))

slots.primerInfo__location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.location, name="primerInfo__location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.primerInfo__location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.microhaplotypesForSample__target_results = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_results, name="microhaplotypesForSample__target_results", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_results'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesForSample__target_results, domain=None, range=Union[Union[dict, MicrohaplotypesForTarget], List[Union[dict, MicrohaplotypesForTarget]]])

slots.microhaplotypeForTarget__reads = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reads, name="microhaplotypeForTarget__reads", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reads'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__reads, domain=None, range=int)

slots.microhaplotypeForTarget__umis = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.umis, name="microhaplotypeForTarget__umis", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('umis'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__umis, domain=None, range=Optional[int])

slots.microhaplotypesForTarget__haps = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.haps, name="microhaplotypesForTarget__haps", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('haps'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesForTarget__haps, domain=None, range=Union[Union[dict, MicrohaplotypeForTarget], List[Union[dict, MicrohaplotypeForTarget]]])

slots.bioinformaticsMethodInfo__demultiplexing_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexing_method, name="bioinformaticsMethodInfo__demultiplexing_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('demultiplexing_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__demultiplexing_method, domain=None, range=Union[dict, BioMethod])

slots.bioinformaticsMethodInfo__denoising_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.denoising_method, name="bioinformaticsMethodInfo__denoising_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('denoising_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__denoising_method, domain=None, range=Union[dict, BioMethod])

slots.bioinformaticsMethodInfo__additional_methods = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.additional_methods, name="bioinformaticsMethodInfo__additional_methods", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('additional_methods'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__additional_methods, domain=None, range=Optional[Union[Union[dict, BioMethod], List[Union[dict, BioMethod]]]])

slots.bioinformaticsMethodInfo__bioinformatics_method_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_method_name, name="bioinformaticsMethodInfo__bioinformatics_method_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_method_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__bioinformatics_method_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.bioMethod__program = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program, name="bioMethod__program", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__program, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.bioMethod__program_description = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_description, name="bioMethod__program_description", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program_description'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__program_description, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.bioMethod__additional_argument = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.additional_argument, name="bioMethod__additional_argument", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('additional_argument'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__additional_argument, domain=None, range=Optional[Union[str, List[str]]],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.experimentInfo__accession = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.accession, name="experimentInfo__accession", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('accession'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experimentInfo__accession, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.experimentInfo__experiment_sample_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_sample_name, name="experimentInfo__experiment_sample_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_sample_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experimentInfo__experiment_sample_name, domain=None, range=URIRef,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__sequencing_info_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_name, name="sequencingInfo__sequencing_info_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__sequencing_info_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__seq_instrument = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_instrument, name="sequencingInfo__seq_instrument", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_instrument'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_instrument, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__seq_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_date, name="sequencingInfo__seq_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_date, domain=None, range=str,
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?'))

slots.sequencingInfo__nucl_acid_ext = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_ext, name="sequencingInfo__nucl_acid_ext", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_ext'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_ext, domain=None, range=Optional[str],
                   pattern=re.compile(r'^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$'))

slots.sequencingInfo__nucl_acid_amp = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_amp, name="sequencingInfo__nucl_acid_amp", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_amp'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_amp, domain=None, range=Optional[str],
                   pattern=re.compile(r'^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$'))

slots.sequencingInfo__nucl_acid_ext_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_ext_date, name="sequencingInfo__nucl_acid_ext_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_ext_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_ext_date, domain=None, range=Optional[str],
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?'))

slots.sequencingInfo__nucl_acid_amp_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_amp_date, name="sequencingInfo__nucl_acid_amp_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_amp_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_amp_date, domain=None, range=Optional[str],
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?'))

slots.sequencingInfo__pcr_cond = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pcr_cond, name="sequencingInfo__pcr_cond", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pcr_cond'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__pcr_cond, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.sequencingInfo__lib_screen = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lib_screen, name="sequencingInfo__lib_screen", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lib_screen'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__lib_screen, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.sequencingInfo__lib_layout = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lib_layout, name="sequencingInfo__lib_layout", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lib_layout'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__lib_layout, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__lib_kit = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lib_kit, name="sequencingInfo__lib_kit", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lib_kit'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__lib_kit, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.sequencingInfo__seq_center = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_center, name="sequencingInfo__seq_center", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_center'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_center, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.parasiteDensity__method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.method, name="parasiteDensity__method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasiteDensity__method, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.parasiteDensity__density = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.density, name="parasiteDensity__density", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('density'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasiteDensity__density, domain=None, range=float,
                   pattern=re.compile(r'^[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?$'))

slots.specimenInfo__specimen_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_name, name="specimenInfo__specimen_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_name, domain=None, range=URIRef,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__samp_taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.samp_taxon_id, name="specimenInfo__samp_taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('samp_taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__samp_taxon_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__individual_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.individual_id, name="specimenInfo__individual_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('individual_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__individual_id, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__host_taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.host_taxon_id, name="specimenInfo__host_taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('host_taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__host_taxon_id, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__alternate_identifiers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alternate_identifiers, name="specimenInfo__alternate_identifiers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alternate_identifiers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__alternate_identifiers, domain=None, range=Optional[Union[str, List[str]]],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__specimen_sex = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_sex, name="specimenInfo__specimen_sex", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_sex'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_sex, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__parasite_density_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasite_density_info, name="specimenInfo__parasite_density_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('parasite_density_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__parasite_density_info, domain=None, range=Optional[Union[Union[dict, ParasiteDensity], List[Union[dict, ParasiteDensity]]]])

slots.specimenInfo__collection_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collection_date, name="specimenInfo__collection_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collection_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collection_date, domain=None, range=str,
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?|NA'))

slots.specimenInfo__date_of_birth = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.date_of_birth, name="specimenInfo__date_of_birth", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('date_of_birth'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__date_of_birth, domain=None, range=Optional[str],
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?|NA'))

slots.specimenInfo__collection_country = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collection_country, name="specimenInfo__collection_country", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collection_country'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collection_country, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__geo_admin1 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin1, name="specimenInfo__geo_admin1", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin1'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin1, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__geo_admin2 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin2, name="specimenInfo__geo_admin2", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin2'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin2, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__geo_admin3 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin3, name="specimenInfo__geo_admin3", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin3'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin3, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__lat_lon = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lat_lon, name="specimenInfo__lat_lon", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lat_lon'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__lat_lon, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[-+]?\d{1,2}(?:\.\d+)?,[-+]?\d{1,3}(?:\.\d+)?$'))

slots.specimenInfo__collector = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collector, name="specimenInfo__collector", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collector'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collector, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__samp_store_loc = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.samp_store_loc, name="specimenInfo__samp_store_loc", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('samp_store_loc'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__samp_store_loc, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__samp_collect_device = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.samp_collect_device, name="specimenInfo__samp_collect_device", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('samp_collect_device'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__samp_collect_device, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__project_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_name, name="specimenInfo__project_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__project_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__sample_comments = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sample_comments, name="specimenInfo__sample_comments", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sample_comments'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__sample_comments, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.bioinformaticsRunInfo__run_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.run_date, name="bioinformaticsRunInfo__run_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('run_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsRunInfo__run_date, domain=None, range=str,
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?'))

slots.bioinformaticsRunInfo__bioinformatics_run_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_run_name, name="bioinformaticsRunInfo__bioinformatics_run_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_run_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsRunInfo__bioinformatics_run_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.pmoGenerationMethod__program_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_name, name="pmoGenerationMethod__program_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pmoGenerationMethod__program_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.pmoHeader__pmo_version = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pmo_version, name="pmoHeader__pmo_version", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pmo_version'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pmoHeader__pmo_version, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.pmoHeader__creation_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.creation_date, name="pmoHeader__creation_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('creation_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pmoHeader__creation_date, domain=None, range=Optional[str],
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?'))

slots.pmoHeader__generation_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.generation_method, name="pmoHeader__generation_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('generation_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pmoHeader__generation_method, domain=None, range=Optional[Union[dict, PmoGenerationMethod]])

slots.stageReadCounts__read_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_count, name="stageReadCounts__read_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.stageReadCounts__read_count, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.stageReadCounts__stage = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.stage, name="stageReadCounts__stage", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('stage'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.stageReadCounts__stage, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.readCountsByStageForTarget__stages = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.stages, name="readCountsByStageForTarget__stages", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('stages'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStageForTarget__stages, domain=None, range=Union[Union[dict, StageReadCounts], List[Union[dict, StageReadCounts]]])

slots.readCountsByStageForExperimentalSample__total_raw_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.total_raw_count, name="readCountsByStageForExperimentalSample__total_raw_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('total_raw_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStageForExperimentalSample__total_raw_count, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.readCountsByStageForExperimentalSample__read_counts_for_targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_counts_for_targets, name="readCountsByStageForExperimentalSample__read_counts_for_targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_counts_for_targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStageForExperimentalSample__read_counts_for_targets, domain=None, range=Optional[Union[Union[dict, ReadCountsByStageForTarget], List[Union[dict, ReadCountsByStageForTarget]]]])

slots.readCountsByStage__read_counts_by_experimental_sample_by_stage = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_counts_by_experimental_sample_by_stage, name="readCountsByStage__read_counts_by_experimental_sample_by_stage", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_counts_by_experimental_sample_by_stage'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStage__read_counts_by_experimental_sample_by_stage, domain=None, range=Union[Union[dict, ReadCountsByStageForExperimentalSample], List[Union[dict, ReadCountsByStageForExperimentalSample]]])

slots.portableMicrohaplotypeObject__experiment_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_info, name="portableMicrohaplotypeObject__experiment_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__experiment_info, domain=None, range=Union[Dict[Union[str, ExperimentInfoExperimentSampleName], Union[dict, ExperimentInfo]], List[Union[dict, ExperimentInfo]]])

slots.portableMicrohaplotypeObject__specimen_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_info, name="portableMicrohaplotypeObject__specimen_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__specimen_info, domain=None, range=Union[Dict[Union[str, SpecimenInfoSpecimenName], Union[dict, SpecimenInfo]], List[Union[dict, SpecimenInfo]]])

slots.portableMicrohaplotypeObject__sequencing_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info, name="portableMicrohaplotypeObject__sequencing_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__sequencing_info, domain=None, range=Union[Union[dict, SequencingInfo], List[Union[dict, SequencingInfo]]])

slots.portableMicrohaplotypeObject__panel_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_info, name="portableMicrohaplotypeObject__panel_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__panel_info, domain=None, range=Union[Union[dict, PanelInfo], List[Union[dict, PanelInfo]]])

slots.portableMicrohaplotypeObject__target_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_info, name="portableMicrohaplotypeObject__target_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__target_info, domain=None, range=Union[Union[dict, TargetInfo], List[Union[dict, TargetInfo]]])

slots.portableMicrohaplotypeObject__targeted_genomes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targeted_genomes, name="portableMicrohaplotypeObject__targeted_genomes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('targeted_genomes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__targeted_genomes, domain=None, range=Union[Union[dict, GenomeInfo], List[Union[dict, GenomeInfo]]])

slots.portableMicrohaplotypeObject__microhaplotypes_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes_info, name="portableMicrohaplotypeObject__microhaplotypes_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__microhaplotypes_info, domain=None, range=Union[dict, RepresentativeMicrohaplotypes])

slots.portableMicrohaplotypeObject__bioinformatics_methods_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_methods_info, name="portableMicrohaplotypeObject__bioinformatics_methods_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_methods_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__bioinformatics_methods_info, domain=None, range=Union[Union[dict, BioinformaticsMethodInfo], List[Union[dict, BioinformaticsMethodInfo]]])

slots.portableMicrohaplotypeObject__bioinformatics_run_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_run_info, name="portableMicrohaplotypeObject__bioinformatics_run_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_run_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__bioinformatics_run_info, domain=None, range=Union[Union[dict, BioinformaticsRunInfo], List[Union[dict, BioinformaticsRunInfo]]])

slots.portableMicrohaplotypeObject__microhaplotypes_detected = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes_detected, name="portableMicrohaplotypeObject__microhaplotypes_detected", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes_detected'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__microhaplotypes_detected, domain=None, range=Union[Union[dict, MicrohaplotypesDetected], List[Union[dict, MicrohaplotypesDetected]]])

slots.portableMicrohaplotypeObject__pmo_header = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pmo_header, name="portableMicrohaplotypeObject__pmo_header", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pmo_header'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__pmo_header, domain=None, range=Union[dict, PmoHeader])

slots.portableMicrohaplotypeObject__read_counts_by_stage = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_counts_by_stage, name="portableMicrohaplotypeObject__read_counts_by_stage", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_counts_by_stage'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__read_counts_by_stage, domain=None, range=Optional[Union[Union[dict, ReadCountsByStage], List[Union[dict, ReadCountsByStage]]]])

slots.RepresentativeMicrohaplotypesForTarget_target_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_id, name="RepresentativeMicrohaplotypesForTarget_target_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypesForTarget_target_id, domain=RepresentativeMicrohaplotypesForTarget, range=Union[int, RepresentativeMicrohaplotypesForTargetTargetId],
                   pattern=re.compile(r'^[0-9]+$'))

slots.ExperimentInfo_specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="ExperimentInfo_specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_specimen_id, domain=ExperimentInfo, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.ExperimentInfo_plate_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, name="ExperimentInfo_plate_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_plate_name, domain=ExperimentInfo, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.ExperimentInfo_plate_row = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, name="ExperimentInfo_plate_row", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_row'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_plate_row, domain=ExperimentInfo, range=Optional[str],
                   pattern=re.compile(r'^[A-z]$'))

slots.ExperimentInfo_plate_col = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, name="ExperimentInfo_plate_col", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_col'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_plate_col, domain=ExperimentInfo, range=Optional[int],
                   pattern=re.compile(r'^[0-9]+$'))