# Auto generated from portable_microhaplotype_object.yaml by pythongen.py version: 0.0.1
# Generation date: 2024-07-22T10:27:24
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
class TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId(extended_str):
    pass


class ExperimentInfoExperimentSampleId(extended_str):
    pass


class SequencingInfoSequencingInfoId(extended_str):
    pass


class SpecimenInfoSpecimenId(extended_str):
    pass


class PortableMicrohaplotypeObjectAnalysisName(extended_str):
    pass


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

    target_id: str = None
    forward_primers: Union[Union[dict, "PrimerInfo"], List[Union[dict, "PrimerInfo"]]] = None
    reverse_primers: Union[Union[dict, "PrimerInfo"], List[Union[dict, "PrimerInfo"]]] = None
    gene_id: Optional[str] = None
    insert_location: Optional[Union[dict, "GenomicLocation"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, str):
            self.target_id = str(self.target_id)

        if self._is_empty(self.forward_primers):
            self.MissingRequiredField("forward_primers")
        self._normalize_inlined_as_dict(slot_name="forward_primers", slot_type=PrimerInfo, key_name="seq", keyed=False)

        if self._is_empty(self.reverse_primers):
            self.MissingRequiredField("reverse_primers")
        self._normalize_inlined_as_dict(slot_name="reverse_primers", slot_type=PrimerInfo, key_name="seq", keyed=False)

        if self.gene_id is not None and not isinstance(self.gene_id, str):
            self.gene_id = str(self.gene_id)

        if self.insert_location is not None and not isinstance(self.insert_location, GenomicLocation):
            self.insert_location = GenomicLocation(**as_dict(self.insert_location))

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

    panel_id: str = None
    target_genome: Union[dict, "GenomeInfo"] = None
    targets: Union[Union[dict, TargetInfo], List[Union[dict, TargetInfo]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.panel_id):
            self.MissingRequiredField("panel_id")
        if not isinstance(self.panel_id, str):
            self.panel_id = str(self.panel_id)

        if self._is_empty(self.target_genome):
            self.MissingRequiredField("target_genome")
        if not isinstance(self.target_genome, GenomeInfo):
            self.target_genome = GenomeInfo(**as_dict(self.target_genome))

        if self._is_empty(self.targets):
            self.MissingRequiredField("targets")
        self._normalize_inlined_as_dict(slot_name="targets", slot_type=TargetInfo, key_name="target_id", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class RepresentativeMicrohaplotypeSequence(YAMLRoot):
    """
    the representative sequence for a microhaplotype, similar to a fast(a/q) format
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotypeSequence"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotypeSequence"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotypeSequence"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypeSequence

    seq: str = None
    microhaplotype_id: str = None
    quality: Optional[str] = None
    pseudocigar: Optional[str] = None
    masking: Optional[Union[Union[dict, "MaskingInfo"], List[Union[dict, "MaskingInfo"]]]] = empty_list()
    alt_annotations: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.seq):
            self.MissingRequiredField("seq")
        if not isinstance(self.seq, str):
            self.seq = str(self.seq)

        if self._is_empty(self.microhaplotype_id):
            self.MissingRequiredField("microhaplotype_id")
        if not isinstance(self.microhaplotype_id, str):
            self.microhaplotype_id = str(self.microhaplotype_id)

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
class RepresentativeMicrohaplotypeSequences(YAMLRoot):
    """
    a list of the representative sequence for a microhaplotypes, similar to a fast(a/q) format
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotypeSequences"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotypeSequences"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotypeSequences"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypeSequences

    target_id: str = None
    seqs: Union[Union[dict, RepresentativeMicrohaplotypeSequence], List[Union[dict, RepresentativeMicrohaplotypeSequence]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, str):
            self.target_id = str(self.target_id)

        if self._is_empty(self.seqs):
            self.MissingRequiredField("seqs")
        if not isinstance(self.seqs, list):
            self.seqs = [self.seqs] if self.seqs is not None else []
        self.seqs = [v if isinstance(v, RepresentativeMicrohaplotypeSequence) else RepresentativeMicrohaplotypeSequence(**as_dict(v)) for v in self.seqs]

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

    tar_amp_bioinformatics_info_id: str = None
    experiment_samples: Union[Union[dict, "MicrohaplotypesForSample"], List[Union[dict, "MicrohaplotypesForSample"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.tar_amp_bioinformatics_info_id):
            self.MissingRequiredField("tar_amp_bioinformatics_info_id")
        if not isinstance(self.tar_amp_bioinformatics_info_id, str):
            self.tar_amp_bioinformatics_info_id = str(self.tar_amp_bioinformatics_info_id)

        if self._is_empty(self.experiment_samples):
            self.MissingRequiredField("experiment_samples")
        self._normalize_inlined_as_dict(slot_name="experiment_samples", slot_type=MicrohaplotypesForSample, key_name="experiment_sample_id", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class DemultiplexedExperimentSamples(YAMLRoot):
    """
    a list of raw reads counts for each experiment sample for all targets within panel
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DemultiplexedExperimentSamples"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DemultiplexedExperimentSamples"
    class_name: ClassVar[str] = "DemultiplexedExperimentSamples"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DemultiplexedExperimentSamples

    tar_amp_bioinformatics_info_id: str = None
    demultiplexed_experiment_samples: Union[Union[dict, "DemultiplexedTargetsForExperimentSample"], List[Union[dict, "DemultiplexedTargetsForExperimentSample"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.tar_amp_bioinformatics_info_id):
            self.MissingRequiredField("tar_amp_bioinformatics_info_id")
        if not isinstance(self.tar_amp_bioinformatics_info_id, str):
            self.tar_amp_bioinformatics_info_id = str(self.tar_amp_bioinformatics_info_id)

        if self._is_empty(self.demultiplexed_experiment_samples):
            self.MissingRequiredField("demultiplexed_experiment_samples")
        self._normalize_inlined_as_dict(slot_name="demultiplexed_experiment_samples", slot_type=DemultiplexedTargetsForExperimentSample, key_name="experiment_sample_id", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class DemultiplexedTargetsForExperimentSample(YAMLRoot):
    """
    a list of raw reads for a experiment sample for all targets within panel
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DemultiplexedTargetsForExperimentSample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DemultiplexedTargetsForExperimentSample"
    class_name: ClassVar[str] = "DemultiplexedTargetsForExperimentSample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DemultiplexedTargetsForExperimentSample

    experiment_sample_id: str = None
    demultiplexed_targets: Union[Union[dict, "DemultiplexedTargetForExperimentSample"], List[Union[dict, "DemultiplexedTargetForExperimentSample"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_sample_id):
            self.MissingRequiredField("experiment_sample_id")
        if not isinstance(self.experiment_sample_id, str):
            self.experiment_sample_id = str(self.experiment_sample_id)

        if self._is_empty(self.demultiplexed_targets):
            self.MissingRequiredField("demultiplexed_targets")
        self._normalize_inlined_as_dict(slot_name="demultiplexed_targets", slot_type=DemultiplexedTargetForExperimentSample, key_name="target_id", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class DemultiplexedTargetForExperimentSample(YAMLRoot):
    """
    the raw read count for a experiment sample for a target
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DemultiplexedTargetForExperimentSample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DemultiplexedTargetForExperimentSample"
    class_name: ClassVar[str] = "DemultiplexedTargetForExperimentSample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DemultiplexedTargetForExperimentSample

    target_id: str = None
    raw_read_count: float = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, str):
            self.target_id = str(self.target_id)

        if self._is_empty(self.raw_read_count):
            self.MissingRequiredField("raw_read_count")
        if not isinstance(self.raw_read_count, float):
            self.raw_read_count = float(self.raw_read_count)

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
    version: str = None
    taxon_id: int = None
    url: str = None
    gff_url: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.name):
            self.MissingRequiredField("name")
        if not isinstance(self.name, str):
            self.name = str(self.name)

        if self._is_empty(self.version):
            self.MissingRequiredField("version")
        if not isinstance(self.version, str):
            self.version = str(self.version)

        if self._is_empty(self.taxon_id):
            self.MissingRequiredField("taxon_id")
        if not isinstance(self.taxon_id, int):
            self.taxon_id = int(self.taxon_id)

        if self._is_empty(self.url):
            self.MissingRequiredField("url")
        if not isinstance(self.url, str):
            self.url = str(self.url)

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

    chrom: str = None
    start: int = None
    end: int = None
    strand: Optional[str] = None
    ref_seq: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
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

    experiment_sample_id: str = None
    target_results: Union[Union[dict, "MicrohaplotypesForTarget"], List[Union[dict, "MicrohaplotypesForTarget"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_sample_id):
            self.MissingRequiredField("experiment_sample_id")
        if not isinstance(self.experiment_sample_id, str):
            self.experiment_sample_id = str(self.experiment_sample_id)

        if self._is_empty(self.target_results):
            self.MissingRequiredField("target_results")
        self._normalize_inlined_as_dict(slot_name="target_results", slot_type=MicrohaplotypesForTarget, key_name="target_id", keyed=False)

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

    microhaplotype_id: str = None
    read_count: float = None
    umi_count: Optional[float] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.microhaplotype_id):
            self.MissingRequiredField("microhaplotype_id")
        if not isinstance(self.microhaplotype_id, str):
            self.microhaplotype_id = str(self.microhaplotype_id)

        if self._is_empty(self.read_count):
            self.MissingRequiredField("read_count")
        if not isinstance(self.read_count, float):
            self.read_count = float(self.read_count)

        if self.umi_count is not None and not isinstance(self.umi_count, float):
            self.umi_count = float(self.umi_count)

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

    target_id: str = None
    microhaplotypes: Union[Union[dict, MicrohaplotypeForTarget], List[Union[dict, MicrohaplotypeForTarget]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, str):
            self.target_id = str(self.target_id)

        if self._is_empty(self.microhaplotypes):
            self.MissingRequiredField("microhaplotypes")
        if not isinstance(self.microhaplotypes, list):
            self.microhaplotypes = [self.microhaplotypes] if self.microhaplotypes is not None else []
        self.microhaplotypes = [v if isinstance(v, MicrohaplotypeForTarget) else MicrohaplotypeForTarget(**as_dict(v)) for v in self.microhaplotypes]

        super().__post_init__(**kwargs)


@dataclass
class TarAmpBioinformaticsInfo(YAMLRoot):
    """
    the targeted amplicon bioinformatics pipeline
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["TarAmpBioinformaticsInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:TarAmpBioinformaticsInfo"
    class_name: ClassVar[str] = "TarAmpBioinformaticsInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.TarAmpBioinformaticsInfo

    tar_amp_bioinformatics_info_id: Union[str, TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId] = None
    demultiplexing_method: Union[dict, "BioMethod"] = None
    denoising_method: Union[dict, "BioMethod"] = None
    population_clustering_method: Union[dict, "BioMethod"] = None
    additional_methods: Optional[Union[Union[dict, "BioMethod"], List[Union[dict, "BioMethod"]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.tar_amp_bioinformatics_info_id):
            self.MissingRequiredField("tar_amp_bioinformatics_info_id")
        if not isinstance(self.tar_amp_bioinformatics_info_id, TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId):
            self.tar_amp_bioinformatics_info_id = TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId(self.tar_amp_bioinformatics_info_id)

        if self._is_empty(self.demultiplexing_method):
            self.MissingRequiredField("demultiplexing_method")
        if not isinstance(self.demultiplexing_method, BioMethod):
            self.demultiplexing_method = BioMethod(**as_dict(self.demultiplexing_method))

        if self._is_empty(self.denoising_method):
            self.MissingRequiredField("denoising_method")
        if not isinstance(self.denoising_method, BioMethod):
            self.denoising_method = BioMethod(**as_dict(self.denoising_method))

        if self._is_empty(self.population_clustering_method):
            self.MissingRequiredField("population_clustering_method")
        if not isinstance(self.population_clustering_method, BioMethod):
            self.population_clustering_method = BioMethod(**as_dict(self.population_clustering_method))

        if not isinstance(self.additional_methods, list):
            self.additional_methods = [self.additional_methods] if self.additional_methods is not None else []
        self.additional_methods = [v if isinstance(v, BioMethod) else BioMethod(**as_dict(v)) for v in self.additional_methods]

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

    program: str = None
    purpose: str = None
    program_version: str = None
    additional_argument: Optional[Union[str, List[str]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.program):
            self.MissingRequiredField("program")
        if not isinstance(self.program, str):
            self.program = str(self.program)

        if self._is_empty(self.purpose):
            self.MissingRequiredField("purpose")
        if not isinstance(self.purpose, str):
            self.purpose = str(self.purpose)

        if self._is_empty(self.program_version):
            self.MissingRequiredField("program_version")
        if not isinstance(self.program_version, str):
            self.program_version = str(self.program_version)

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

    experiment_sample_id: Union[str, ExperimentInfoExperimentSampleId] = None
    sequencing_info_id: str = None
    specimen_id: str = None
    panel_id: str = None
    plate_name: Optional[str] = None
    plate_row: Optional[str] = None
    plate_col: Optional[int] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_sample_id):
            self.MissingRequiredField("experiment_sample_id")
        if not isinstance(self.experiment_sample_id, ExperimentInfoExperimentSampleId):
            self.experiment_sample_id = ExperimentInfoExperimentSampleId(self.experiment_sample_id)

        if self._is_empty(self.sequencing_info_id):
            self.MissingRequiredField("sequencing_info_id")
        if not isinstance(self.sequencing_info_id, str):
            self.sequencing_info_id = str(self.sequencing_info_id)

        if self._is_empty(self.specimen_id):
            self.MissingRequiredField("specimen_id")
        if not isinstance(self.specimen_id, str):
            self.specimen_id = str(self.specimen_id)

        if self._is_empty(self.panel_id):
            self.MissingRequiredField("panel_id")
        if not isinstance(self.panel_id, str):
            self.panel_id = str(self.panel_id)

        if self.plate_name is not None and not isinstance(self.plate_name, str):
            self.plate_name = str(self.plate_name)

        if self.plate_row is not None and not isinstance(self.plate_row, str):
            self.plate_row = str(self.plate_row)

        if self.plate_col is not None and not isinstance(self.plate_col, int):
            self.plate_col = int(self.plate_col)

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

    sequencing_info_id: Union[str, SequencingInfoSequencingInfoId] = None
    seq_instrument: str = None
    seq_date: str = None
    nucl_acid_ext: str = None
    nucl_acid_amp: str = None
    nucl_acid_ext_date: str = None
    nucl_acid_amp_date: str = None
    pcr_cond: str = None
    lib_screen: str = None
    lib_layout: str = None
    lib_kit: str = None
    seq_center: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.sequencing_info_id):
            self.MissingRequiredField("sequencing_info_id")
        if not isinstance(self.sequencing_info_id, SequencingInfoSequencingInfoId):
            self.sequencing_info_id = SequencingInfoSequencingInfoId(self.sequencing_info_id)

        if self._is_empty(self.seq_instrument):
            self.MissingRequiredField("seq_instrument")
        if not isinstance(self.seq_instrument, str):
            self.seq_instrument = str(self.seq_instrument)

        if self._is_empty(self.seq_date):
            self.MissingRequiredField("seq_date")
        if not isinstance(self.seq_date, str):
            self.seq_date = str(self.seq_date)

        if self._is_empty(self.nucl_acid_ext):
            self.MissingRequiredField("nucl_acid_ext")
        if not isinstance(self.nucl_acid_ext, str):
            self.nucl_acid_ext = str(self.nucl_acid_ext)

        if self._is_empty(self.nucl_acid_amp):
            self.MissingRequiredField("nucl_acid_amp")
        if not isinstance(self.nucl_acid_amp, str):
            self.nucl_acid_amp = str(self.nucl_acid_amp)

        if self._is_empty(self.nucl_acid_ext_date):
            self.MissingRequiredField("nucl_acid_ext_date")
        if not isinstance(self.nucl_acid_ext_date, str):
            self.nucl_acid_ext_date = str(self.nucl_acid_ext_date)

        if self._is_empty(self.nucl_acid_amp_date):
            self.MissingRequiredField("nucl_acid_amp_date")
        if not isinstance(self.nucl_acid_amp_date, str):
            self.nucl_acid_amp_date = str(self.nucl_acid_amp_date)

        if self._is_empty(self.pcr_cond):
            self.MissingRequiredField("pcr_cond")
        if not isinstance(self.pcr_cond, str):
            self.pcr_cond = str(self.pcr_cond)

        if self._is_empty(self.lib_screen):
            self.MissingRequiredField("lib_screen")
        if not isinstance(self.lib_screen, str):
            self.lib_screen = str(self.lib_screen)

        if self._is_empty(self.lib_layout):
            self.MissingRequiredField("lib_layout")
        if not isinstance(self.lib_layout, str):
            self.lib_layout = str(self.lib_layout)

        if self._is_empty(self.lib_kit):
            self.MissingRequiredField("lib_kit")
        if not isinstance(self.lib_kit, str):
            self.lib_kit = str(self.lib_kit)

        if self._is_empty(self.seq_center):
            self.MissingRequiredField("seq_center")
        if not isinstance(self.seq_center, str):
            self.seq_center = str(self.seq_center)

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

    specimen_id: Union[str, SpecimenInfoSpecimenId] = None
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
    individual_id: Optional[str] = None
    host_taxon_id: Optional[int] = None
    alternate_identifiers: Optional[Union[str, List[str]]] = empty_list()
    parasite_density: Optional[int] = None
    geo_admin1: Optional[str] = None
    geo_admin2: Optional[str] = None
    geo_admin3: Optional[str] = None
    lat_lon: Optional[str] = None
    accession: Optional[str] = None
    sample_comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.specimen_id):
            self.MissingRequiredField("specimen_id")
        if not isinstance(self.specimen_id, SpecimenInfoSpecimenId):
            self.specimen_id = SpecimenInfoSpecimenId(self.specimen_id)

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

        if self.individual_id is not None and not isinstance(self.individual_id, str):
            self.individual_id = str(self.individual_id)

        if self.host_taxon_id is not None and not isinstance(self.host_taxon_id, int):
            self.host_taxon_id = int(self.host_taxon_id)

        if not isinstance(self.alternate_identifiers, list):
            self.alternate_identifiers = [self.alternate_identifiers] if self.alternate_identifiers is not None else []
        self.alternate_identifiers = [v if isinstance(v, str) else str(v) for v in self.alternate_identifiers]

        if self.parasite_density is not None and not isinstance(self.parasite_density, int):
            self.parasite_density = int(self.parasite_density)

        if self.geo_admin1 is not None and not isinstance(self.geo_admin1, str):
            self.geo_admin1 = str(self.geo_admin1)

        if self.geo_admin2 is not None and not isinstance(self.geo_admin2, str):
            self.geo_admin2 = str(self.geo_admin2)

        if self.geo_admin3 is not None and not isinstance(self.geo_admin3, str):
            self.geo_admin3 = str(self.geo_admin3)

        if self.lat_lon is not None and not isinstance(self.lat_lon, str):
            self.lat_lon = str(self.lat_lon)

        if self.accession is not None and not isinstance(self.accession, str):
            self.accession = str(self.accession)

        if self.sample_comments is not None and not isinstance(self.sample_comments, str):
            self.sample_comments = str(self.sample_comments)

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

    analysis_name: Union[str, PortableMicrohaplotypeObjectAnalysisName] = None
    experiment_infos: Union[Dict[Union[str, ExperimentInfoExperimentSampleId], Union[dict, ExperimentInfo]], List[Union[dict, ExperimentInfo]]] = empty_dict()
    specimen_infos: Union[Dict[Union[str, SpecimenInfoSpecimenId], Union[dict, SpecimenInfo]], List[Union[dict, SpecimenInfo]]] = empty_dict()
    sequencing_infos: Union[dict, SequencingInfo] = None
    panel_info: Union[dict, PanelInfo] = None
    representative_microhaplotype_sequences: Union[Union[dict, RepresentativeMicrohaplotypeSequences], List[Union[dict, RepresentativeMicrohaplotypeSequences]]] = None
    microhaplotypes_detected: Union[dict, MicrohaplotypesDetected] = None
    taramp_bioinformatics_infos: Union[dict, TarAmpBioinformaticsInfo] = None
    target_demultiplexed_experiment_samples: Optional[Union[dict, DemultiplexedExperimentSamples]] = None
    postprocessing_bioinformatics_infos: Optional[Union[dict, BioMethod]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.analysis_name):
            self.MissingRequiredField("analysis_name")
        if not isinstance(self.analysis_name, PortableMicrohaplotypeObjectAnalysisName):
            self.analysis_name = PortableMicrohaplotypeObjectAnalysisName(self.analysis_name)

        if self._is_empty(self.experiment_infos):
            self.MissingRequiredField("experiment_infos")
        self._normalize_inlined_as_list(slot_name="experiment_infos", slot_type=ExperimentInfo, key_name="experiment_sample_id", keyed=True)

        if self._is_empty(self.specimen_infos):
            self.MissingRequiredField("specimen_infos")
        self._normalize_inlined_as_list(slot_name="specimen_infos", slot_type=SpecimenInfo, key_name="specimen_id", keyed=True)

        if self._is_empty(self.sequencing_infos):
            self.MissingRequiredField("sequencing_infos")
        if not isinstance(self.sequencing_infos, SequencingInfo):
            self.sequencing_infos = SequencingInfo(**as_dict(self.sequencing_infos))

        if self._is_empty(self.panel_info):
            self.MissingRequiredField("panel_info")
        if not isinstance(self.panel_info, PanelInfo):
            self.panel_info = PanelInfo(**as_dict(self.panel_info))

        if self._is_empty(self.representative_microhaplotype_sequences):
            self.MissingRequiredField("representative_microhaplotype_sequences")
        self._normalize_inlined_as_dict(slot_name="representative_microhaplotype_sequences", slot_type=RepresentativeMicrohaplotypeSequences, key_name="target_id", keyed=False)

        if self._is_empty(self.microhaplotypes_detected):
            self.MissingRequiredField("microhaplotypes_detected")
        if not isinstance(self.microhaplotypes_detected, MicrohaplotypesDetected):
            self.microhaplotypes_detected = MicrohaplotypesDetected(**as_dict(self.microhaplotypes_detected))

        if self._is_empty(self.taramp_bioinformatics_infos):
            self.MissingRequiredField("taramp_bioinformatics_infos")
        if not isinstance(self.taramp_bioinformatics_infos, TarAmpBioinformaticsInfo):
            self.taramp_bioinformatics_infos = TarAmpBioinformaticsInfo(**as_dict(self.taramp_bioinformatics_infos))

        if self.target_demultiplexed_experiment_samples is not None and not isinstance(self.target_demultiplexed_experiment_samples, DemultiplexedExperimentSamples):
            self.target_demultiplexed_experiment_samples = DemultiplexedExperimentSamples(**as_dict(self.target_demultiplexed_experiment_samples))

        if self.postprocessing_bioinformatics_infos is not None and not isinstance(self.postprocessing_bioinformatics_infos, BioMethod):
            self.postprocessing_bioinformatics_infos = BioMethod(**as_dict(self.postprocessing_bioinformatics_infos))

        super().__post_init__(**kwargs)


# Enumerations


# Slots
class slots:
    pass

slots.seq = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq, name="seq", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq, domain=None, range=str,
                   pattern=re.compile(r'^[A-z]$'))

slots.specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.experiment_sample_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_sample_id, name="experiment_sample_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_sample_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_sample_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.target_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_id, name="target_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.microhaplotype_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotype_id, name="microhaplotype_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotype_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotype_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.panel_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_id, name="panel_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.plate_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, name="plate_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.plate_row = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, name="plate_row", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_row'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z]$'))

slots.plate_col = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, name="plate_col", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_col'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]$'))

slots.sequencing_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, name="sequencing_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.tar_amp_bioinformatics_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tar_amp_bioinformatics_info_id, name="tar_amp_bioinformatics_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('tar_amp_bioinformatics_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tar_amp_bioinformatics_info_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.targetInfo__gene_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gene_id, name="targetInfo__gene_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gene_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__gene_id, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.targetInfo__insert_location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.insert_location, name="targetInfo__insert_location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('insert_location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__insert_location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.targetInfo__forward_primers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.forward_primers, name="targetInfo__forward_primers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('forward_primers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__forward_primers, domain=None, range=Union[Union[dict, PrimerInfo], List[Union[dict, PrimerInfo]]])

slots.targetInfo__reverse_primers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reverse_primers, name="targetInfo__reverse_primers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reverse_primers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__reverse_primers, domain=None, range=Union[Union[dict, PrimerInfo], List[Union[dict, PrimerInfo]]])

slots.panelInfo__target_genome = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_genome, name="panelInfo__target_genome", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_genome'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__target_genome, domain=None, range=Union[dict, GenomeInfo])

slots.panelInfo__targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targets, name="panelInfo__targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__targets, domain=None, range=Union[Union[dict, TargetInfo], List[Union[dict, TargetInfo]]])

slots.representativeMicrohaplotypeSequence__quality = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.quality, name="representativeMicrohaplotypeSequence__quality", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('quality'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequence__quality, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.representativeMicrohaplotypeSequence__pseudocigar = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar, name="representativeMicrohaplotypeSequence__pseudocigar", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pseudocigar'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequence__pseudocigar, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.representativeMicrohaplotypeSequence__masking = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.masking, name="representativeMicrohaplotypeSequence__masking", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('masking'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequence__masking, domain=None, range=Optional[Union[Union[dict, MaskingInfo], List[Union[dict, MaskingInfo]]]])

slots.representativeMicrohaplotypeSequence__alt_annotations = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alt_annotations, name="representativeMicrohaplotypeSequence__alt_annotations", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alt_annotations'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequence__alt_annotations, domain=None, range=Optional[Union[str, List[str]]])

slots.maskingInfo__seq_start = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_start, name="maskingInfo__seq_start", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_start'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__seq_start, domain=None, range=int)

slots.maskingInfo__seq_segment_size = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_segment_size, name="maskingInfo__seq_segment_size", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_segment_size'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__seq_segment_size, domain=None, range=int)

slots.maskingInfo__replacement_size = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.replacement_size, name="maskingInfo__replacement_size", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('replacement_size'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__replacement_size, domain=None, range=int)

slots.representativeMicrohaplotypeSequences__seqs = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seqs, name="representativeMicrohaplotypeSequences__seqs", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seqs'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequences__seqs, domain=None, range=Union[Union[dict, RepresentativeMicrohaplotypeSequence], List[Union[dict, RepresentativeMicrohaplotypeSequence]]])

slots.microhaplotypesDetected__experiment_samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_samples, name="microhaplotypesDetected__experiment_samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesDetected__experiment_samples, domain=None, range=Union[Union[dict, MicrohaplotypesForSample], List[Union[dict, MicrohaplotypesForSample]]])

slots.demultiplexedExperimentSamples__demultiplexed_experiment_samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexed_experiment_samples, name="demultiplexedExperimentSamples__demultiplexed_experiment_samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('demultiplexed_experiment_samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexedExperimentSamples__demultiplexed_experiment_samples, domain=None, range=Union[Union[dict, DemultiplexedTargetsForExperimentSample], List[Union[dict, DemultiplexedTargetsForExperimentSample]]])

slots.demultiplexedTargetsForExperimentSample__demultiplexed_targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexed_targets, name="demultiplexedTargetsForExperimentSample__demultiplexed_targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('demultiplexed_targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexedTargetsForExperimentSample__demultiplexed_targets, domain=None, range=Union[Union[dict, DemultiplexedTargetForExperimentSample], List[Union[dict, DemultiplexedTargetForExperimentSample]]])

slots.demultiplexedTargetForExperimentSample__raw_read_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.raw_read_count, name="demultiplexedTargetForExperimentSample__raw_read_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('raw_read_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexedTargetForExperimentSample__raw_read_count, domain=None, range=float)

slots.genomeInfo__name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.name, name="genomeInfo__name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.genomeInfo__version = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.version, name="genomeInfo__version", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('version'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__version, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.genomeInfo__taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.taxon_id, name="genomeInfo__taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__taxon_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]$'))

slots.genomeInfo__url = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.url, name="genomeInfo__url", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('url'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__url, domain=None, range=str,
                   pattern=re.compile(r'r"^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$"'))

slots.genomeInfo__gff_url = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gff_url, name="genomeInfo__gff_url", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gff_url'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__gff_url, domain=None, range=Optional[str],
                   pattern=re.compile(r'r"^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$"'))

slots.genomicLocation__chrom = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.chrom, name="genomicLocation__chrom", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('chrom'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__chrom, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.genomicLocation__start = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.start, name="genomicLocation__start", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('start'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__start, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]$'))

slots.genomicLocation__end = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.end, name="genomicLocation__end", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('end'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__end, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]$'))

slots.genomicLocation__strand = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.strand, name="genomicLocation__strand", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('strand'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__strand, domain=None, range=Optional[str],
                   pattern=re.compile(r'r'[+-]''))

slots.genomicLocation__ref_seq = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ref_seq, name="genomicLocation__ref_seq", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('ref_seq'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__ref_seq, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-]$'))

slots.primerInfo__location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.location, name="primerInfo__location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.primerInfo__location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.microhaplotypesForSample__target_results = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_results, name="microhaplotypesForSample__target_results", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_results'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesForSample__target_results, domain=None, range=Union[Union[dict, MicrohaplotypesForTarget], List[Union[dict, MicrohaplotypesForTarget]]])

slots.microhaplotypeForTarget__read_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_count, name="microhaplotypeForTarget__read_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__read_count, domain=None, range=float)

slots.microhaplotypeForTarget__umi_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.umi_count, name="microhaplotypeForTarget__umi_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('umi_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__umi_count, domain=None, range=Optional[float])

slots.microhaplotypesForTarget__microhaplotypes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes, name="microhaplotypesForTarget__microhaplotypes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesForTarget__microhaplotypes, domain=None, range=Union[Union[dict, MicrohaplotypeForTarget], List[Union[dict, MicrohaplotypeForTarget]]])

slots.tarAmpBioinformaticsInfo__demultiplexing_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexing_method, name="tarAmpBioinformaticsInfo__demultiplexing_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('demultiplexing_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tarAmpBioinformaticsInfo__demultiplexing_method, domain=None, range=Union[dict, BioMethod])

slots.tarAmpBioinformaticsInfo__denoising_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.denoising_method, name="tarAmpBioinformaticsInfo__denoising_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('denoising_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tarAmpBioinformaticsInfo__denoising_method, domain=None, range=Union[dict, BioMethod])

slots.tarAmpBioinformaticsInfo__population_clustering_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.population_clustering_method, name="tarAmpBioinformaticsInfo__population_clustering_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('population_clustering_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tarAmpBioinformaticsInfo__population_clustering_method, domain=None, range=Union[dict, BioMethod])

slots.tarAmpBioinformaticsInfo__additional_methods = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.additional_methods, name="tarAmpBioinformaticsInfo__additional_methods", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('additional_methods'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tarAmpBioinformaticsInfo__additional_methods, domain=None, range=Optional[Union[Union[dict, BioMethod], List[Union[dict, BioMethod]]]])

slots.bioMethod__program = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program, name="bioMethod__program", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__program, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.bioMethod__purpose = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.purpose, name="bioMethod__purpose", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('purpose'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__purpose, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.bioMethod__program_version = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_version, name="bioMethod__program_version", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program_version'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__program_version, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.bioMethod__additional_argument = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.additional_argument, name="bioMethod__additional_argument", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('additional_argument'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__additional_argument, domain=None, range=Optional[Union[str, List[str]]],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.sequencingInfo__seq_instrument = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_instrument, name="sequencingInfo__seq_instrument", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_instrument'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_instrument, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.sequencingInfo__seq_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_date, name="sequencingInfo__seq_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_date, domain=None, range=str,
                   pattern=re.compile(r'r"\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?"'))

slots.sequencingInfo__nucl_acid_ext = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_ext, name="sequencingInfo__nucl_acid_ext", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_ext'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_ext, domain=None, range=str,
                   pattern=re.compile(r'r"^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$"'))

slots.sequencingInfo__nucl_acid_amp = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_amp, name="sequencingInfo__nucl_acid_amp", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_amp'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_amp, domain=None, range=str,
                   pattern=re.compile(r'r"^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$"'))

slots.sequencingInfo__nucl_acid_ext_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_ext_date, name="sequencingInfo__nucl_acid_ext_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_ext_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_ext_date, domain=None, range=str,
                   pattern=re.compile(r'r"\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?"'))

slots.sequencingInfo__nucl_acid_amp_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.nucl_acid_amp_date, name="sequencingInfo__nucl_acid_amp_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('nucl_acid_amp_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__nucl_acid_amp_date, domain=None, range=str,
                   pattern=re.compile(r'r"\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?"'))

slots.sequencingInfo__pcr_cond = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pcr_cond, name="sequencingInfo__pcr_cond", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pcr_cond'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__pcr_cond, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.sequencingInfo__lib_screen = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lib_screen, name="sequencingInfo__lib_screen", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lib_screen'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__lib_screen, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.sequencingInfo__lib_layout = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lib_layout, name="sequencingInfo__lib_layout", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lib_layout'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__lib_layout, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.sequencingInfo__lib_kit = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lib_kit, name="sequencingInfo__lib_kit", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lib_kit'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__lib_kit, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.sequencingInfo__seq_center = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_center, name="sequencingInfo__seq_center", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_center'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_center, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__samp_taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.samp_taxon_id, name="specimenInfo__samp_taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('samp_taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__samp_taxon_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]$'))

slots.specimenInfo__individual_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.individual_id, name="specimenInfo__individual_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('individual_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__individual_id, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__host_taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.host_taxon_id, name="specimenInfo__host_taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('host_taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__host_taxon_id, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]$'))

slots.specimenInfo__alternate_identifiers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alternate_identifiers, name="specimenInfo__alternate_identifiers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alternate_identifiers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__alternate_identifiers, domain=None, range=Optional[Union[str, List[str]]],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__parasite_density = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasite_density, name="specimenInfo__parasite_density", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('parasite_density'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__parasite_density, domain=None, range=Optional[int],
                   pattern=re.compile(r'r'^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$''))

slots.specimenInfo__collection_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collection_date, name="specimenInfo__collection_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collection_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collection_date, domain=None, range=str,
                   pattern=re.compile(r'r"\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?"'))

slots.specimenInfo__collection_country = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collection_country, name="specimenInfo__collection_country", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collection_country'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collection_country, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__geo_admin1 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin1, name="specimenInfo__geo_admin1", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin1'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin1, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__geo_admin2 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin2, name="specimenInfo__geo_admin2", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin2'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin2, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__geo_admin3 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin3, name="specimenInfo__geo_admin3", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin3'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin3, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__lat_lon = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lat_lon, name="specimenInfo__lat_lon", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lat_lon'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__lat_lon, domain=None, range=Optional[str],
                   pattern=re.compile(r'r'^[-+]?\d{1,2}(?:\.\d+)?,[-+]?\d{1,3}(?:\.\d+)?$''))

slots.specimenInfo__collector = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collector, name="specimenInfo__collector", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collector'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collector, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__samp_store_loc = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.samp_store_loc, name="specimenInfo__samp_store_loc", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('samp_store_loc'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__samp_store_loc, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__samp_collect_device = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.samp_collect_device, name="specimenInfo__samp_collect_device", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('samp_collect_device'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__samp_collect_device, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__project_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_name, name="specimenInfo__project_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__project_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.specimenInfo__accession = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.accession, name="specimenInfo__accession", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('accession'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__accession, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.specimenInfo__sample_comments = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sample_comments, name="specimenInfo__sample_comments", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sample_comments'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__sample_comments, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.portableMicrohaplotypeObject__analysis_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.analysis_name, name="portableMicrohaplotypeObject__analysis_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('analysis_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__analysis_name, domain=None, range=URIRef)

slots.portableMicrohaplotypeObject__experiment_infos = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_infos, name="portableMicrohaplotypeObject__experiment_infos", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_infos'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__experiment_infos, domain=None, range=Union[Dict[Union[str, ExperimentInfoExperimentSampleId], Union[dict, ExperimentInfo]], List[Union[dict, ExperimentInfo]]])

slots.portableMicrohaplotypeObject__specimen_infos = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_infos, name="portableMicrohaplotypeObject__specimen_infos", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_infos'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__specimen_infos, domain=None, range=Union[Dict[Union[str, SpecimenInfoSpecimenId], Union[dict, SpecimenInfo]], List[Union[dict, SpecimenInfo]]])

slots.portableMicrohaplotypeObject__sequencing_infos = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_infos, name="portableMicrohaplotypeObject__sequencing_infos", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_infos'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__sequencing_infos, domain=None, range=Union[dict, SequencingInfo])

slots.portableMicrohaplotypeObject__panel_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_info, name="portableMicrohaplotypeObject__panel_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__panel_info, domain=None, range=Union[dict, PanelInfo])

slots.portableMicrohaplotypeObject__representative_microhaplotype_sequences = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representative_microhaplotype_sequences, name="portableMicrohaplotypeObject__representative_microhaplotype_sequences", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('representative_microhaplotype_sequences'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__representative_microhaplotype_sequences, domain=None, range=Union[Union[dict, RepresentativeMicrohaplotypeSequences], List[Union[dict, RepresentativeMicrohaplotypeSequences]]])

slots.portableMicrohaplotypeObject__microhaplotypes_detected = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes_detected, name="portableMicrohaplotypeObject__microhaplotypes_detected", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes_detected'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__microhaplotypes_detected, domain=None, range=Union[dict, MicrohaplotypesDetected])

slots.portableMicrohaplotypeObject__target_demultiplexed_experiment_samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_demultiplexed_experiment_samples, name="portableMicrohaplotypeObject__target_demultiplexed_experiment_samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_demultiplexed_experiment_samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__target_demultiplexed_experiment_samples, domain=None, range=Optional[Union[dict, DemultiplexedExperimentSamples]])

slots.portableMicrohaplotypeObject__taramp_bioinformatics_infos = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.taramp_bioinformatics_infos, name="portableMicrohaplotypeObject__taramp_bioinformatics_infos", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('taramp_bioinformatics_infos'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__taramp_bioinformatics_infos, domain=None, range=Union[dict, TarAmpBioinformaticsInfo])

slots.portableMicrohaplotypeObject__postprocessing_bioinformatics_infos = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.postprocessing_bioinformatics_infos, name="portableMicrohaplotypeObject__postprocessing_bioinformatics_infos", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('postprocessing_bioinformatics_infos'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__postprocessing_bioinformatics_infos, domain=None, range=Optional[Union[dict, BioMethod]])

slots.TarAmpBioinformaticsInfo_tar_amp_bioinformatics_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tar_amp_bioinformatics_info_id, name="TarAmpBioinformaticsInfo_tar_amp_bioinformatics_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('tar_amp_bioinformatics_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.TarAmpBioinformaticsInfo_tar_amp_bioinformatics_info_id, domain=TarAmpBioinformaticsInfo, range=Union[str, TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.ExperimentInfo_specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="ExperimentInfo_specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_specimen_id, domain=ExperimentInfo, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.ExperimentInfo_panel_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_id, name="ExperimentInfo_panel_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_panel_id, domain=ExperimentInfo, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.ExperimentInfo_experiment_sample_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_sample_id, name="ExperimentInfo_experiment_sample_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_sample_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_experiment_sample_id, domain=ExperimentInfo, range=Union[str, ExperimentInfoExperimentSampleId],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.ExperimentInfo_plate_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, name="ExperimentInfo_plate_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_plate_name, domain=ExperimentInfo, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.ExperimentInfo_plate_row = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, name="ExperimentInfo_plate_row", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_row'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_plate_row, domain=ExperimentInfo, range=Optional[str],
                   pattern=re.compile(r'^[A-z]$'))

slots.ExperimentInfo_plate_col = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, name="ExperimentInfo_plate_col", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_col'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_plate_col, domain=ExperimentInfo, range=Optional[int],
                   pattern=re.compile(r'^[0-9]$'))

slots.SequencingInfo_sequencing_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, name="SequencingInfo_sequencing_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.SequencingInfo_sequencing_info_id, domain=SequencingInfo, range=Union[str, SequencingInfoSequencingInfoId],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.SpecimenInfo_specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="SpecimenInfo_specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.SpecimenInfo_specimen_id, domain=SpecimenInfo, range=Union[str, SpecimenInfoSpecimenId],
                   pattern=re.compile(r'^[A-z-._0-9]$'))