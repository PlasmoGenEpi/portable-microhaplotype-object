# Auto generated from portable_microhaplotype_object.yaml by pythongen.py version: 0.0.1
# Generation date: 2024-07-15T17:04:13
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


class ExperimentInfoExperimentId(extended_str):
    pass


class SequencingInfoSequencingInfoId(extended_str):
    pass


class SpecimenInfoSpecimenId(extended_str):
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
    insert_location: Union[dict, "GenomicLocation"] = None
    forward_primers: Union[dict, "Primers"] = None
    reverse_primers: Union[dict, "Primers"] = None
    gene_id: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, str):
            self.target_id = str(self.target_id)

        if self._is_empty(self.insert_location):
            self.MissingRequiredField("insert_location")
        if not isinstance(self.insert_location, GenomicLocation):
            self.insert_location = GenomicLocation(**as_dict(self.insert_location))

        if self._is_empty(self.forward_primers):
            self.MissingRequiredField("forward_primers")
        if not isinstance(self.forward_primers, Primers):
            self.forward_primers = Primers(**as_dict(self.forward_primers))

        if self._is_empty(self.reverse_primers):
            self.MissingRequiredField("reverse_primers")
        if not isinstance(self.reverse_primers, Primers):
            self.reverse_primers = Primers(**as_dict(self.reverse_primers))

        if self.gene_id is not None and not isinstance(self.gene_id, str):
            self.gene_id = str(self.gene_id)

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

        if not isinstance(self.alt_annotations, list):
            self.alt_annotations = [self.alt_annotations] if self.alt_annotations is not None else []
        self.alt_annotations = [v if isinstance(v, str) else str(v) for v in self.alt_annotations]

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

    sequencing_id: str = None
    bioinformatics_id: str = None
    samples: Union[Union[dict, "MicrohaplotypesForSample"], List[Union[dict, "MicrohaplotypesForSample"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.sequencing_id):
            self.MissingRequiredField("sequencing_id")
        if not isinstance(self.sequencing_id, str):
            self.sequencing_id = str(self.sequencing_id)

        if self._is_empty(self.bioinformatics_id):
            self.MissingRequiredField("bioinformatics_id")
        if not isinstance(self.bioinformatics_id, str):
            self.bioinformatics_id = str(self.bioinformatics_id)

        if self._is_empty(self.samples):
            self.MissingRequiredField("samples")
        self._normalize_inlined_as_dict(slot_name="samples", slot_type=MicrohaplotypesForSample, key_name="experiment_id", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class DemultiplexedSamples(YAMLRoot):
    """
    a list of raw reads counts for each sample for all targets within panel
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DemultiplexedSamples"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DemultiplexedSamples"
    class_name: ClassVar[str] = "DemultiplexedSamples"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DemultiplexedSamples

    demultiplexed_samples: Union[Union[dict, "DemultiplexedTargetsForSample"], List[Union[dict, "DemultiplexedTargetsForSample"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.demultiplexed_samples):
            self.MissingRequiredField("demultiplexed_samples")
        self._normalize_inlined_as_dict(slot_name="demultiplexed_samples", slot_type=DemultiplexedTargetsForSample, key_name="experiment_id", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class DemultiplexedTargetsForSample(YAMLRoot):
    """
    a list of raw reads for a sample for all targets within panel
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DemultiplexedTargetsForSample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DemultiplexedTargetsForSample"
    class_name: ClassVar[str] = "DemultiplexedTargetsForSample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DemultiplexedTargetsForSample

    experiment_id: str = None
    demultiplexed_targets: Union[Union[dict, "DemultiplexedTargetForSample"], List[Union[dict, "DemultiplexedTargetForSample"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_id):
            self.MissingRequiredField("experiment_id")
        if not isinstance(self.experiment_id, str):
            self.experiment_id = str(self.experiment_id)

        if self._is_empty(self.demultiplexed_targets):
            self.MissingRequiredField("demultiplexed_targets")
        self._normalize_inlined_as_dict(slot_name="demultiplexed_targets", slot_type=DemultiplexedTargetForSample, key_name="target_id", keyed=False)

        super().__post_init__(**kwargs)


@dataclass
class DemultiplexedTargetForSample(YAMLRoot):
    """
    the raw read count for a sample for a target
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DemultiplexedTargetForSample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DemultiplexedTargetForSample"
    class_name: ClassVar[str] = "DemultiplexedTargetForSample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DemultiplexedTargetForSample

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
    location: Union[dict, GenomicLocation] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.seq):
            self.MissingRequiredField("seq")
        if not isinstance(self.seq, str):
            self.seq = str(self.seq)

        if self._is_empty(self.location):
            self.MissingRequiredField("location")
        if not isinstance(self.location, GenomicLocation):
            self.location = GenomicLocation(**as_dict(self.location))

        super().__post_init__(**kwargs)


@dataclass
class Primers(YAMLRoot):
    """
    A holder of primer sequences
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["Primers"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:Primers"
    class_name: ClassVar[str] = "Primers"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.Primers

    entries: Optional[Union[Union[dict, PrimerInfo], List[Union[dict, PrimerInfo]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        self._normalize_inlined_as_dict(slot_name="entries", slot_type=PrimerInfo, key_name="seq", keyed=False)

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

    experiment_id: str = None
    target_results: Union[Union[dict, "MicrohaplotypesForTarget"], List[Union[dict, "MicrohaplotypesForTarget"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_id):
            self.MissingRequiredField("experiment_id")
        if not isinstance(self.experiment_id, str):
            self.experiment_id = str(self.experiment_id)

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

    experiment_id: Union[str, ExperimentInfoExperimentId] = None
    sequencing_info_id: str = None
    plate_name: str = None
    plate_row: str = None
    plate_col: str = None
    specimen_id: str = None
    panel_id: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_id):
            self.MissingRequiredField("experiment_id")
        if not isinstance(self.experiment_id, ExperimentInfoExperimentId):
            self.experiment_id = ExperimentInfoExperimentId(self.experiment_id)

        if self._is_empty(self.sequencing_info_id):
            self.MissingRequiredField("sequencing_info_id")
        if not isinstance(self.sequencing_info_id, str):
            self.sequencing_info_id = str(self.sequencing_info_id)

        if self._is_empty(self.plate_name):
            self.MissingRequiredField("plate_name")
        if not isinstance(self.plate_name, str):
            self.plate_name = str(self.plate_name)

        if self._is_empty(self.plate_row):
            self.MissingRequiredField("plate_row")
        if not isinstance(self.plate_row, str):
            self.plate_row = str(self.plate_row)

        if self._is_empty(self.plate_col):
            self.MissingRequiredField("plate_col")
        if not isinstance(self.plate_col, str):
            self.plate_col = str(self.plate_col)

        if self._is_empty(self.specimen_id):
            self.MissingRequiredField("specimen_id")
        if not isinstance(self.specimen_id, str):
            self.specimen_id = str(self.specimen_id)

        if self._is_empty(self.panel_id):
            self.MissingRequiredField("panel_id")
        if not isinstance(self.panel_id, str):
            self.panel_id = str(self.panel_id)

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
    plate_name: str = None
    plate_row: str = None
    plate_col: str = None
    samp_taxon_id: int = None
    collection_date: str = None
    lat_lon: str = None
    geo_loc_name: str = None
    collector: str = None
    samp_store_loc: str = None
    samp_collect_device: str = None
    project_name: str = None
    individual_id: Optional[str] = None
    host_taxon_id: Optional[int] = None
    alternate_identifiers: Optional[Union[str, List[str]]] = empty_list()
    parasite_density: Optional[int] = None
    accession: Optional[str] = None
    sample_comments: Optional[str] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.specimen_id):
            self.MissingRequiredField("specimen_id")
        if not isinstance(self.specimen_id, SpecimenInfoSpecimenId):
            self.specimen_id = SpecimenInfoSpecimenId(self.specimen_id)

        if self._is_empty(self.plate_name):
            self.MissingRequiredField("plate_name")
        if not isinstance(self.plate_name, str):
            self.plate_name = str(self.plate_name)

        if self._is_empty(self.plate_row):
            self.MissingRequiredField("plate_row")
        if not isinstance(self.plate_row, str):
            self.plate_row = str(self.plate_row)

        if self._is_empty(self.plate_col):
            self.MissingRequiredField("plate_col")
        if not isinstance(self.plate_col, str):
            self.plate_col = str(self.plate_col)

        if self._is_empty(self.samp_taxon_id):
            self.MissingRequiredField("samp_taxon_id")
        if not isinstance(self.samp_taxon_id, int):
            self.samp_taxon_id = int(self.samp_taxon_id)

        if self._is_empty(self.collection_date):
            self.MissingRequiredField("collection_date")
        if not isinstance(self.collection_date, str):
            self.collection_date = str(self.collection_date)

        if self._is_empty(self.lat_lon):
            self.MissingRequiredField("lat_lon")
        if not isinstance(self.lat_lon, str):
            self.lat_lon = str(self.lat_lon)

        if self._is_empty(self.geo_loc_name):
            self.MissingRequiredField("geo_loc_name")
        if not isinstance(self.geo_loc_name, str):
            self.geo_loc_name = str(self.geo_loc_name)

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

        if self.individual_id is not None and not isinstance(self.individual_id, str):
            self.individual_id = str(self.individual_id)

        if self.host_taxon_id is not None and not isinstance(self.host_taxon_id, int):
            self.host_taxon_id = int(self.host_taxon_id)

        if not isinstance(self.alternate_identifiers, list):
            self.alternate_identifiers = [self.alternate_identifiers] if self.alternate_identifiers is not None else []
        self.alternate_identifiers = [v if isinstance(v, str) else str(v) for v in self.alternate_identifiers]

        if self.parasite_density is not None and not isinstance(self.parasite_density, int):
            self.parasite_density = int(self.parasite_density)

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

    experiment_infos: Union[Dict[Union[str, ExperimentInfoExperimentId], Union[dict, ExperimentInfo]], List[Union[dict, ExperimentInfo]]] = empty_dict()
    specimen_infos: Union[Dict[Union[str, SpecimenInfoSpecimenId], Union[dict, SpecimenInfo]], List[Union[dict, SpecimenInfo]]] = empty_dict()
    sequencing_info: Union[dict, SequencingInfo] = None
    bioinformatics_info: Union[str, TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId] = None
    panel_info: Union[dict, PanelInfo] = None
    representative_microhaplotype_sequences: Union[Union[dict, RepresentativeMicrohaplotypeSequences], List[Union[dict, RepresentativeMicrohaplotypeSequences]]] = None
    microhaplotypes_detected: Union[dict, MicrohaplotypesDetected] = None
    target_demultiplexed_samples: Union[dict, DemultiplexedSamples] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.experiment_infos):
            self.MissingRequiredField("experiment_infos")
        self._normalize_inlined_as_list(slot_name="experiment_infos", slot_type=ExperimentInfo, key_name="experiment_id", keyed=True)

        if self._is_empty(self.specimen_infos):
            self.MissingRequiredField("specimen_infos")
        self._normalize_inlined_as_list(slot_name="specimen_infos", slot_type=SpecimenInfo, key_name="specimen_id", keyed=True)

        if self._is_empty(self.sequencing_info):
            self.MissingRequiredField("sequencing_info")
        if not isinstance(self.sequencing_info, SequencingInfo):
            self.sequencing_info = SequencingInfo(**as_dict(self.sequencing_info))

        if self._is_empty(self.bioinformatics_info):
            self.MissingRequiredField("bioinformatics_info")
        if not isinstance(self.bioinformatics_info, TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId):
            self.bioinformatics_info = TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId(self.bioinformatics_info)

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

        if self._is_empty(self.target_demultiplexed_samples):
            self.MissingRequiredField("target_demultiplexed_samples")
        if not isinstance(self.target_demultiplexed_samples, DemultiplexedSamples):
            self.target_demultiplexed_samples = DemultiplexedSamples(**as_dict(self.target_demultiplexed_samples))

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

slots.experiment_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_id, name="experiment_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_id, domain=None, range=str,
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
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.plate_row = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, name="plate_row", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_row'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.plate_col = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, name="plate_col", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_col'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.sequencing_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, name="sequencing_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.targetInfo__gene_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gene_id, name="targetInfo__gene_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gene_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__gene_id, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.targetInfo__insert_location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.insert_location, name="targetInfo__insert_location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('insert_location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__insert_location, domain=None, range=Union[dict, GenomicLocation])

slots.targetInfo__forward_primers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.forward_primers, name="targetInfo__forward_primers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('forward_primers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__forward_primers, domain=None, range=Union[dict, Primers])

slots.targetInfo__reverse_primers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reverse_primers, name="targetInfo__reverse_primers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reverse_primers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__reverse_primers, domain=None, range=Union[dict, Primers])

slots.panelInfo__target_genome = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_genome, name="panelInfo__target_genome", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_genome'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__target_genome, domain=None, range=Union[dict, GenomeInfo])

slots.panelInfo__targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targets, name="panelInfo__targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__targets, domain=None, range=Union[Union[dict, TargetInfo], List[Union[dict, TargetInfo]]])

slots.representativeMicrohaplotypeSequence__quality = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.quality, name="representativeMicrohaplotypeSequence__quality", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('quality'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequence__quality, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.representativeMicrohaplotypeSequence__alt_annotations = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alt_annotations, name="representativeMicrohaplotypeSequence__alt_annotations", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alt_annotations'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequence__alt_annotations, domain=None, range=Optional[Union[str, List[str]]])

slots.representativeMicrohaplotypeSequences__seqs = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seqs, name="representativeMicrohaplotypeSequences__seqs", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seqs'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypeSequences__seqs, domain=None, range=Union[Union[dict, RepresentativeMicrohaplotypeSequence], List[Union[dict, RepresentativeMicrohaplotypeSequence]]])

slots.microhaplotypesDetected__sequencing_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_id, name="microhaplotypesDetected__sequencing_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesDetected__sequencing_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.microhaplotypesDetected__bioinformatics_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_id, name="microhaplotypesDetected__bioinformatics_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesDetected__bioinformatics_id, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.microhaplotypesDetected__samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.samples, name="microhaplotypesDetected__samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesDetected__samples, domain=None, range=Union[Union[dict, MicrohaplotypesForSample], List[Union[dict, MicrohaplotypesForSample]]])

slots.demultiplexedSamples__demultiplexed_samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexed_samples, name="demultiplexedSamples__demultiplexed_samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('demultiplexed_samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexedSamples__demultiplexed_samples, domain=None, range=Union[Union[dict, DemultiplexedTargetsForSample], List[Union[dict, DemultiplexedTargetsForSample]]])

slots.demultiplexedTargetsForSample__demultiplexed_targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexed_targets, name="demultiplexedTargetsForSample__demultiplexed_targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('demultiplexed_targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexedTargetsForSample__demultiplexed_targets, domain=None, range=Union[Union[dict, DemultiplexedTargetForSample], List[Union[dict, DemultiplexedTargetForSample]]])

slots.demultiplexedTargetForSample__raw_read_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.raw_read_count, name="demultiplexedTargetForSample__raw_read_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('raw_read_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexedTargetForSample__raw_read_count, domain=None, range=float)

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

slots.primerInfo__location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.location, name="primerInfo__location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.primerInfo__location, domain=None, range=Union[dict, GenomicLocation])

slots.primers__entries = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.entries, name="primers__entries", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('entries'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.primers__entries, domain=None, range=Optional[Union[Union[dict, PrimerInfo], List[Union[dict, PrimerInfo]]]])

slots.microhaplotypesForSample__target_results = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_results, name="microhaplotypesForSample__target_results", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_results'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesForSample__target_results, domain=None, range=Union[Union[dict, MicrohaplotypesForTarget], List[Union[dict, MicrohaplotypesForTarget]]])

slots.microhaplotypeForTarget__read_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_count, name="microhaplotypeForTarget__read_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__read_count, domain=None, range=float)

slots.microhaplotypeForTarget__umi_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.umi_count, name="microhaplotypeForTarget__umi_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('umi_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__umi_count, domain=None, range=Optional[float])

slots.microhaplotypesForTarget__microhaplotypes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes, name="microhaplotypesForTarget__microhaplotypes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypesForTarget__microhaplotypes, domain=None, range=Union[Union[dict, MicrohaplotypeForTarget], List[Union[dict, MicrohaplotypeForTarget]]])

slots.tarAmpBioinformaticsInfo__tar_amp_bioinformatics_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tar_amp_bioinformatics_info_id, name="tarAmpBioinformaticsInfo__tar_amp_bioinformatics_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('tar_amp_bioinformatics_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.tarAmpBioinformaticsInfo__tar_amp_bioinformatics_info_id, domain=None, range=URIRef,
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

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

slots.specimenInfo__lat_lon = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lat_lon, name="specimenInfo__lat_lon", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lat_lon'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__lat_lon, domain=None, range=str,
                   pattern=re.compile(r'r'^[-+]?\d{1,2}(?:\.\d+)?,[-+]?\d{1,3}(?:\.\d+)?$''))

slots.specimenInfo__geo_loc_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_loc_name, name="specimenInfo__geo_loc_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_loc_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_loc_name, domain=None, range=str,
                   pattern=re.compile(r'^([^\s-]{1,2}|[^\s-]+.+[^\s-]+): ([^\s-]{1,2}|[^\s-]+.+[^\s-]+), ([^\s-]{1,2}|[^\s-]+.+[^\s-]+)$'))

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

slots.portableMicrohaplotypeObject__experiment_infos = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_infos, name="portableMicrohaplotypeObject__experiment_infos", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_infos'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__experiment_infos, domain=None, range=Union[Dict[Union[str, ExperimentInfoExperimentId], Union[dict, ExperimentInfo]], List[Union[dict, ExperimentInfo]]])

slots.portableMicrohaplotypeObject__specimen_infos = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_infos, name="portableMicrohaplotypeObject__specimen_infos", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_infos'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__specimen_infos, domain=None, range=Union[Dict[Union[str, SpecimenInfoSpecimenId], Union[dict, SpecimenInfo]], List[Union[dict, SpecimenInfo]]])

slots.portableMicrohaplotypeObject__sequencing_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info, name="portableMicrohaplotypeObject__sequencing_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__sequencing_info, domain=None, range=Union[dict, SequencingInfo])

slots.portableMicrohaplotypeObject__bioinformatics_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_info, name="portableMicrohaplotypeObject__bioinformatics_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__bioinformatics_info, domain=None, range=Union[str, TarAmpBioinformaticsInfoTarAmpBioinformaticsInfoId])

slots.portableMicrohaplotypeObject__panel_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_info, name="portableMicrohaplotypeObject__panel_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__panel_info, domain=None, range=Union[dict, PanelInfo])

slots.portableMicrohaplotypeObject__representative_microhaplotype_sequences = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representative_microhaplotype_sequences, name="portableMicrohaplotypeObject__representative_microhaplotype_sequences", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('representative_microhaplotype_sequences'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__representative_microhaplotype_sequences, domain=None, range=Union[Union[dict, RepresentativeMicrohaplotypeSequences], List[Union[dict, RepresentativeMicrohaplotypeSequences]]])

slots.portableMicrohaplotypeObject__microhaplotypes_detected = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes_detected, name="portableMicrohaplotypeObject__microhaplotypes_detected", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes_detected'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__microhaplotypes_detected, domain=None, range=Union[dict, MicrohaplotypesDetected])

slots.portableMicrohaplotypeObject__target_demultiplexed_samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_demultiplexed_samples, name="portableMicrohaplotypeObject__target_demultiplexed_samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_demultiplexed_samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__target_demultiplexed_samples, domain=None, range=Union[dict, DemultiplexedSamples])

slots.ExperimentInfo_specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="ExperimentInfo_specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_specimen_id, domain=ExperimentInfo, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.ExperimentInfo_panel_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_id, name="ExperimentInfo_panel_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_panel_id, domain=ExperimentInfo, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]$'))

slots.ExperimentInfo_experiment_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.experiment_id, name="ExperimentInfo_experiment_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('experiment_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ExperimentInfo_experiment_id, domain=ExperimentInfo, range=Union[str, ExperimentInfoExperimentId],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.SequencingInfo_sequencing_info_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_id, name="SequencingInfo_sequencing_info_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.SequencingInfo_sequencing_info_id, domain=SequencingInfo, range=Union[str, SequencingInfoSequencingInfoId],
                   pattern=re.compile(r'^[A-z-._0-9 ]$'))

slots.SpecimenInfo_specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="SpecimenInfo_specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.SpecimenInfo_specimen_id, domain=SpecimenInfo, range=Union[str, SpecimenInfoSpecimenId],
                   pattern=re.compile(r'^[A-z-._0-9]$'))