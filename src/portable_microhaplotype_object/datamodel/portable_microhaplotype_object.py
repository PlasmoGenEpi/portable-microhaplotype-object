# Auto generated from portable_microhaplotype_object.yaml by pythongen.py version: 0.0.1
# Generation date: 2025-07-16T16:07:42
# Schema: portable-microhaplotype-object
#
# id: https://plasmogenepi.github.io/portable-microhaplotype-object
# description: A schema to define the minimum amount of data needed to export a microhaplotype calling pipeline analysis with associated metadata
# license: GNU GPL v3.0

import dataclasses
import re
from dataclasses import dataclass
from datetime import (
    date,
    datetime,
    time
)
from typing import (
    Any,
    ClassVar,
    Dict,
    List,
    Optional,
    Union
)

from jsonasobj2 import (
    JsonObj,
    as_dict
)
from linkml_runtime.linkml_model.meta import (
    EnumDefinition,
    PermissibleValue,
    PvFormulaOptions
)
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from linkml_runtime.utils.formatutils import (
    camelcase,
    sfx,
    underscore
)
from linkml_runtime.utils.metamodelcore import (
    bnode,
    empty_dict,
    empty_list
)
from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.yamlutils import (
    YAMLRoot,
    extended_float,
    extended_int,
    extended_str
)
from rdflib import (
    Namespace,
    URIRef
)

from linkml_runtime.linkml_model.types import Boolean, Double, Integer, String
from linkml_runtime.utils.metamodelcore import Bool

metamodel_version = "1.7.0"
version = None

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


class LibrarySampleInfoLibrarySampleName(extended_str):
    pass


class ProjectInfoProjectName(extended_str):
    pass


class SpecimenInfoSpecimenName(extended_str):
    pass


@dataclass(repr=False)
class MarkerOfInterest(YAMLRoot):
    """
    A specific genomic location of interest, e.g. drug resistance, or other phenotypical marker
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MarkerOfInterest"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MarkerOfInterest"
    class_name: ClassVar[str] = "MarkerOfInterest"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MarkerOfInterest

    marker_location: Union[dict, "GenomicLocation"] = None
    associations: Optional[Union[str, list[str]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.marker_location):
            self.MissingRequiredField("marker_location")
        if not isinstance(self.marker_location, GenomicLocation):
            self.marker_location = GenomicLocation(**as_dict(self.marker_location))

        if not isinstance(self.associations, list):
            self.associations = [self.associations] if self.associations is not None else []
        self.associations = [v if isinstance(v, str) else str(v) for v in self.associations]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class TargetInfo(YAMLRoot):
    """
    Information about a specific target within a genome
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["TargetInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:TargetInfo"
    class_name: ClassVar[str] = "TargetInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.TargetInfo

    target_name: str = None
    forward_primer: Union[dict, "PrimerInfo"] = None
    reverse_primer: Union[dict, "PrimerInfo"] = None
    gene_name: Optional[str] = None
    insert_location: Optional[Union[dict, "GenomicLocation"]] = None
    markers_of_interest: Optional[Union[Union[dict, MarkerOfInterest], list[Union[dict, MarkerOfInterest]]]] = empty_list()
    target_attributes: Optional[Union[str, list[str]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.target_name):
            self.MissingRequiredField("target_name")
        if not isinstance(self.target_name, str):
            self.target_name = str(self.target_name)

        if self._is_empty(self.forward_primer):
            self.MissingRequiredField("forward_primer")
        if not isinstance(self.forward_primer, PrimerInfo):
            self.forward_primer = PrimerInfo(**as_dict(self.forward_primer))

        if self._is_empty(self.reverse_primer):
            self.MissingRequiredField("reverse_primer")
        if not isinstance(self.reverse_primer, PrimerInfo):
            self.reverse_primer = PrimerInfo(**as_dict(self.reverse_primer))

        if self.gene_name is not None and not isinstance(self.gene_name, str):
            self.gene_name = str(self.gene_name)

        if self.insert_location is not None and not isinstance(self.insert_location, GenomicLocation):
            self.insert_location = GenomicLocation(**as_dict(self.insert_location))

        self._normalize_inlined_as_dict(slot_name="markers_of_interest", slot_type=MarkerOfInterest, key_name="marker_location", keyed=False)

        if not isinstance(self.target_attributes, list):
            self.target_attributes = [self.target_attributes] if self.target_attributes is not None else []
        self.target_attributes = [v if isinstance(v, str) else str(v) for v in self.target_attributes]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ReactionInfo(YAMLRoot):
    """
    information on a panel of targeted amplicon primer pairs
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReactionInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReactionInfo"
    class_name: ClassVar[str] = "ReactionInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReactionInfo

    panel_targets: Union[int, list[int]] = None
    reaction_name: str = None

    def __post_init__(self, *_: str, **kwargs: Any):
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


@dataclass(repr=False)
class PanelInfo(YAMLRoot):
    """
    information on a panel of targeted amplicon primer pairs
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PanelInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PanelInfo"
    class_name: ClassVar[str] = "PanelInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PanelInfo

    reactions: Union[Union[dict, ReactionInfo], list[Union[dict, ReactionInfo]]] = None
    panel_name: str = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.reactions):
            self.MissingRequiredField("reactions")
        self._normalize_inlined_as_dict(slot_name="reactions", slot_type=ReactionInfo, key_name="panel_targets", keyed=False)

        if self._is_empty(self.panel_name):
            self.MissingRequiredField("panel_name")
        if not isinstance(self.panel_name, str):
            self.panel_name = str(self.panel_name)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Pseudocigar(YAMLRoot):
    """
    information on pseudocigar for a sequence
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["Pseudocigar"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:Pseudocigar"
    class_name: ClassVar[str] = "Pseudocigar"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.Pseudocigar

    pseudocigar_seq: str = None
    ref_loc: Union[dict, "GenomicLocation"] = None
    pseudocigar_generation_description: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.pseudocigar_seq):
            self.MissingRequiredField("pseudocigar_seq")
        if not isinstance(self.pseudocigar_seq, str):
            self.pseudocigar_seq = str(self.pseudocigar_seq)

        if self._is_empty(self.ref_loc):
            self.MissingRequiredField("ref_loc")
        if not isinstance(self.ref_loc, GenomicLocation):
            self.ref_loc = GenomicLocation(**as_dict(self.ref_loc))

        if self.pseudocigar_generation_description is not None and not isinstance(self.pseudocigar_generation_description, str):
            self.pseudocigar_generation_description = str(self.pseudocigar_generation_description)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class RepresentativeMicrohaplotype(YAMLRoot):
    """
    the representative sequence for a microhaplotype, similar to a fast(a/q) format
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotype"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotype"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotype"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotype

    seq: str = None
    microhaplotype_name: Optional[str] = None
    quality: Optional[str] = None
    pseudocigar: Optional[Union[dict, Pseudocigar]] = None
    masking: Optional[Union[Union[dict, "MaskingInfo"], list[Union[dict, "MaskingInfo"]]]] = empty_list()
    alt_annotations: Optional[Union[str, list[str]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.seq):
            self.MissingRequiredField("seq")
        if not isinstance(self.seq, str):
            self.seq = str(self.seq)

        if self.microhaplotype_name is not None and not isinstance(self.microhaplotype_name, str):
            self.microhaplotype_name = str(self.microhaplotype_name)

        if self.quality is not None and not isinstance(self.quality, str):
            self.quality = str(self.quality)

        if self.pseudocigar is not None and not isinstance(self.pseudocigar, Pseudocigar):
            self.pseudocigar = Pseudocigar(**as_dict(self.pseudocigar))

        self._normalize_inlined_as_dict(slot_name="masking", slot_type=MaskingInfo, key_name="seq_start", keyed=False)

        if not isinstance(self.alt_annotations, list):
            self.alt_annotations = [self.alt_annotations] if self.alt_annotations is not None else []
        self.alt_annotations = [v if isinstance(v, str) else str(v) for v in self.alt_annotations]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class MaskingInfo(YAMLRoot):
    """
    information about a subsegment of the sequence that should be masked
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MaskingInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MaskingInfo"
    class_name: ClassVar[str] = "MaskingInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MaskingInfo

    seq_start: int = None
    seq_segment_size: int = None
    replacement_size: int = None
    masking_generation_description: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
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

        if self.masking_generation_description is not None and not isinstance(self.masking_generation_description, str):
            self.masking_generation_description = str(self.masking_generation_description)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class RepresentativeMicrohaplotypes(YAMLRoot):
    """
    a collection of representative sequences for microhaplotypes for all targets
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotypes"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotypes"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotypes"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypes

    targets: Union[dict[Union[int, RepresentativeMicrohaplotypesForTargetTargetId], Union[dict, "RepresentativeMicrohaplotypesForTarget"]], list[Union[dict, "RepresentativeMicrohaplotypesForTarget"]]] = empty_dict()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.targets):
            self.MissingRequiredField("targets")
        self._normalize_inlined_as_list(slot_name="targets", slot_type=RepresentativeMicrohaplotypesForTarget, key_name="target_id", keyed=True)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class RepresentativeMicrohaplotypesForTarget(YAMLRoot):
    """
    a list of the representative sequence for a microhaplotypes, similar to a fast(a/q) format
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["RepresentativeMicrohaplotypesForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:RepresentativeMicrohaplotypesForTarget"
    class_name: ClassVar[str] = "RepresentativeMicrohaplotypesForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypesForTarget

    target_id: Union[int, RepresentativeMicrohaplotypesForTargetTargetId] = None
    microhaplotypes: Union[Union[dict, RepresentativeMicrohaplotype], list[Union[dict, RepresentativeMicrohaplotype]]] = None
    mhap_location: Optional[Union[dict, "GenomicLocation"]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.target_id):
            self.MissingRequiredField("target_id")
        if not isinstance(self.target_id, RepresentativeMicrohaplotypesForTargetTargetId):
            self.target_id = RepresentativeMicrohaplotypesForTargetTargetId(self.target_id)

        if self._is_empty(self.microhaplotypes):
            self.MissingRequiredField("microhaplotypes")
        self._normalize_inlined_as_dict(slot_name="microhaplotypes", slot_type=RepresentativeMicrohaplotype, key_name="seq", keyed=False)

        if self.mhap_location is not None and not isinstance(self.mhap_location, GenomicLocation):
            self.mhap_location = GenomicLocation(**as_dict(self.mhap_location))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class DetectedMicrohaplotypes(YAMLRoot):
    """
    the microhaplotypes detected in a targeted amplicon analysis
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DetectedMicrohaplotypes"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DetectedMicrohaplotypes"
    class_name: ClassVar[str] = "DetectedMicrohaplotypes"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DetectedMicrohaplotypes

    bioinformatics_run_id: int = None
    library_samples: Union[Union[dict, "DetectedMicrohaplotypesForSample"], list[Union[dict, "DetectedMicrohaplotypesForSample"]]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.bioinformatics_run_id):
            self.MissingRequiredField("bioinformatics_run_id")
        if not isinstance(self.bioinformatics_run_id, int):
            self.bioinformatics_run_id = int(self.bioinformatics_run_id)

        if self._is_empty(self.library_samples):
            self.MissingRequiredField("library_samples")
        if not isinstance(self.library_samples, list):
            self.library_samples = [self.library_samples] if self.library_samples is not None else []
        self.library_samples = [v if isinstance(v, DetectedMicrohaplotypesForSample) else DetectedMicrohaplotypesForSample(**as_dict(v)) for v in self.library_samples]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class GenomeInfo(YAMLRoot):
    """
    information on a genome
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["GenomeInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:GenomeInfo"
    class_name: ClassVar[str] = "GenomeInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.GenomeInfo

    name: str = None
    genome_version: str = None
    taxon_id: Union[int, list[int]] = None
    url: str = None
    chromosomes: Optional[Union[str, list[str]]] = empty_list()
    gff_url: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
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
        if not isinstance(self.taxon_id, list):
            self.taxon_id = [self.taxon_id] if self.taxon_id is not None else []
        self.taxon_id = [v if isinstance(v, int) else int(v) for v in self.taxon_id]

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


@dataclass(repr=False)
class GenomicLocation(YAMLRoot):
    """
    information on the genomic location of specific sequence
    """
    _inherited_slots: ClassVar[list[str]] = []

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
    alt_seq: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
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

        if self.alt_seq is not None and not isinstance(self.alt_seq, str):
            self.alt_seq = str(self.alt_seq)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class PrimerInfo(YAMLRoot):
    """
    information on a primer sequence
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PrimerInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PrimerInfo"
    class_name: ClassVar[str] = "PrimerInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PrimerInfo

    seq: str = None
    location: Optional[Union[dict, GenomicLocation]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.seq):
            self.MissingRequiredField("seq")
        if not isinstance(self.seq, str):
            self.seq = str(self.seq)

        if self.location is not None and not isinstance(self.location, GenomicLocation):
            self.location = GenomicLocation(**as_dict(self.location))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class DetectedMicrohaplotypesForSample(YAMLRoot):
    """
    Microhaplotypes detected for a sample for all targets
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DetectedMicrohaplotypesForSample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DetectedMicrohaplotypesForSample"
    class_name: ClassVar[str] = "DetectedMicrohaplotypesForSample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DetectedMicrohaplotypesForSample

    library_sample_id: int = None
    target_results: Union[Union[dict, "DetectedMicrohaplotypesForTarget"], list[Union[dict, "DetectedMicrohaplotypesForTarget"]]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.library_sample_id):
            self.MissingRequiredField("library_sample_id")
        if not isinstance(self.library_sample_id, int):
            self.library_sample_id = int(self.library_sample_id)

        if self._is_empty(self.target_results):
            self.MissingRequiredField("target_results")
        if not isinstance(self.target_results, list):
            self.target_results = [self.target_results] if self.target_results is not None else []
        self.target_results = [v if isinstance(v, DetectedMicrohaplotypesForTarget) else DetectedMicrohaplotypesForTarget(**as_dict(v)) for v in self.target_results]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class MicrohaplotypeForTarget(YAMLRoot):
    """
    Microhaplotype detected for a specific target
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["MicrohaplotypeForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:MicrohaplotypeForTarget"
    class_name: ClassVar[str] = "MicrohaplotypeForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.MicrohaplotypeForTarget

    mhap_id: int = None
    reads: int = None
    umis: Optional[int] = None

    def __post_init__(self, *_: str, **kwargs: Any):
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


@dataclass(repr=False)
class DetectedMicrohaplotypesForTarget(YAMLRoot):
    """
    Microhaplotypes detected for a specific target
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["DetectedMicrohaplotypesForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:DetectedMicrohaplotypesForTarget"
    class_name: ClassVar[str] = "DetectedMicrohaplotypesForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.DetectedMicrohaplotypesForTarget

    mhaps_target_id: int = None
    mhaps: Union[Union[dict, MicrohaplotypeForTarget], list[Union[dict, MicrohaplotypeForTarget]]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.mhaps_target_id):
            self.MissingRequiredField("mhaps_target_id")
        if not isinstance(self.mhaps_target_id, int):
            self.mhaps_target_id = int(self.mhaps_target_id)

        if self._is_empty(self.mhaps):
            self.MissingRequiredField("mhaps")
        if not isinstance(self.mhaps, list):
            self.mhaps = [self.mhaps] if self.mhaps is not None else []
        self.mhaps = [v if isinstance(v, MicrohaplotypeForTarget) else MicrohaplotypeForTarget(**as_dict(v)) for v in self.mhaps]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class BioinformaticsMethodInfo(YAMLRoot):
    """
    the targeted amplicon bioinformatics pipeline
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["BioinformaticsMethodInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:BioinformaticsMethodInfo"
    class_name: ClassVar[str] = "BioinformaticsMethodInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.BioinformaticsMethodInfo

    demultiplexing_method: Union[dict, "BioMethod"] = None
    denoising_method: Union[dict, "BioMethod"] = None
    additional_methods: Optional[Union[Union[dict, "BioMethod"], list[Union[dict, "BioMethod"]]]] = empty_list()
    bioinformatics_method_name: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
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


@dataclass(repr=False)
class BioMethod(YAMLRoot):
    """
    methodology description of a portion of a bioinformatics pipeline
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["BioMethod"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:BioMethod"
    class_name: ClassVar[str] = "BioMethod"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.BioMethod

    program_version: str = None
    program: str = None
    program_description: Optional[str] = None
    program_url: Optional[str] = None
    additional_argument: Optional[Union[str, list[str]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
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

        if self.program_url is not None and not isinstance(self.program_url, str):
            self.program_url = str(self.program_url)

        if not isinstance(self.additional_argument, list):
            self.additional_argument = [self.additional_argument] if self.additional_argument is not None else []
        self.additional_argument = [v if isinstance(v, str) else str(v) for v in self.additional_argument]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class PlateInfo(YAMLRoot):
    """
    Information about a plate location in a standard 96 well plate
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PlateInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PlateInfo"
    class_name: ClassVar[str] = "PlateInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PlateInfo

    plate_name: Optional[str] = None
    plate_row: Optional[str] = None
    plate_col: Optional[int] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.plate_name is not None and not isinstance(self.plate_name, str):
            self.plate_name = str(self.plate_name)

        if self.plate_row is not None and not isinstance(self.plate_row, str):
            self.plate_row = str(self.plate_row)

        if self.plate_col is not None and not isinstance(self.plate_col, int):
            self.plate_col = int(self.plate_col)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class LibrarySampleInfo(YAMLRoot):
    """
    Information about a specific amplification and sequencing of a specimen
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["LibrarySampleInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:LibrarySampleInfo"
    class_name: ClassVar[str] = "LibrarySampleInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.LibrarySampleInfo

    library_sample_name: Union[str, LibrarySampleInfoLibrarySampleName] = None
    sequencing_info_id: int = None
    specimen_id: int = None
    panel_id: int = None
    fastqs_loc: Optional[str] = None
    run_accession: Optional[str] = None
    library_prep_plate_info: Optional[Union[dict, PlateInfo]] = None
    qpcr_parasite_density_info: Optional[Union[Union[dict, "ParasiteDensity"], list[Union[dict, "ParasiteDensity"]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.library_sample_name):
            self.MissingRequiredField("library_sample_name")
        if not isinstance(self.library_sample_name, LibrarySampleInfoLibrarySampleName):
            self.library_sample_name = LibrarySampleInfoLibrarySampleName(self.library_sample_name)

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

        if self.fastqs_loc is not None and not isinstance(self.fastqs_loc, str):
            self.fastqs_loc = str(self.fastqs_loc)

        if self.run_accession is not None and not isinstance(self.run_accession, str):
            self.run_accession = str(self.run_accession)

        if self.library_prep_plate_info is not None and not isinstance(self.library_prep_plate_info, PlateInfo):
            self.library_prep_plate_info = PlateInfo(**as_dict(self.library_prep_plate_info))

        if not isinstance(self.qpcr_parasite_density_info, list):
            self.qpcr_parasite_density_info = [self.qpcr_parasite_density_info] if self.qpcr_parasite_density_info is not None else []
        self.qpcr_parasite_density_info = [v if isinstance(v, ParasiteDensity) else ParasiteDensity(**as_dict(v)) for v in self.qpcr_parasite_density_info]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class SequencingInfo(YAMLRoot):
    """
    Information on sequencing info
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["SequencingInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:SequencingInfo"
    class_name: ClassVar[str] = "SequencingInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.SequencingInfo

    sequencing_info_name: str = None
    seq_platform: str = None
    seq_instrument_model: str = None
    library_layout: str = None
    library_strategy: str = None
    library_source: str = None
    library_selection: str = None
    seq_date: Optional[str] = None
    nucl_acid_ext: Optional[str] = None
    nucl_acid_amp: Optional[str] = None
    nucl_acid_ext_date: Optional[str] = None
    nucl_acid_amp_date: Optional[str] = None
    pcr_cond: Optional[str] = None
    library_screen: Optional[str] = None
    library_kit: Optional[str] = None
    seq_center: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.sequencing_info_name):
            self.MissingRequiredField("sequencing_info_name")
        if not isinstance(self.sequencing_info_name, str):
            self.sequencing_info_name = str(self.sequencing_info_name)

        if self._is_empty(self.seq_platform):
            self.MissingRequiredField("seq_platform")
        if not isinstance(self.seq_platform, str):
            self.seq_platform = str(self.seq_platform)

        if self._is_empty(self.seq_instrument_model):
            self.MissingRequiredField("seq_instrument_model")
        if not isinstance(self.seq_instrument_model, str):
            self.seq_instrument_model = str(self.seq_instrument_model)

        if self._is_empty(self.library_layout):
            self.MissingRequiredField("library_layout")
        if not isinstance(self.library_layout, str):
            self.library_layout = str(self.library_layout)

        if self._is_empty(self.library_strategy):
            self.MissingRequiredField("library_strategy")
        if not isinstance(self.library_strategy, str):
            self.library_strategy = str(self.library_strategy)

        if self._is_empty(self.library_source):
            self.MissingRequiredField("library_source")
        if not isinstance(self.library_source, str):
            self.library_source = str(self.library_source)

        if self._is_empty(self.library_selection):
            self.MissingRequiredField("library_selection")
        if not isinstance(self.library_selection, str):
            self.library_selection = str(self.library_selection)

        if self.seq_date is not None and not isinstance(self.seq_date, str):
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

        if self.library_screen is not None and not isinstance(self.library_screen, str):
            self.library_screen = str(self.library_screen)

        if self.library_kit is not None and not isinstance(self.library_kit, str):
            self.library_kit = str(self.library_kit)

        if self.seq_center is not None and not isinstance(self.seq_center, str):
            self.seq_center = str(self.seq_center)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ParasiteDensity(YAMLRoot):
    """
    method and value of determined parasite density
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ParasiteDensity"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ParasiteDensity"
    class_name: ClassVar[str] = "ParasiteDensity"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ParasiteDensity

    density_method: str = None
    parasite_density: float = None
    date_measured: Optional[str] = None
    density_method_comments: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.density_method):
            self.MissingRequiredField("density_method")
        if not isinstance(self.density_method, str):
            self.density_method = str(self.density_method)

        if self._is_empty(self.parasite_density):
            self.MissingRequiredField("parasite_density")
        if not isinstance(self.parasite_density, float):
            self.parasite_density = float(self.parasite_density)

        if self.date_measured is not None and not isinstance(self.date_measured, str):
            self.date_measured = str(self.date_measured)

        if self.density_method_comments is not None and not isinstance(self.density_method_comments, str):
            self.density_method_comments = str(self.density_method_comments)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ProjectInfo(YAMLRoot):
    """
    Information on project info
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ProjectInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ProjectInfo"
    class_name: ClassVar[str] = "ProjectInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ProjectInfo

    project_name: Union[str, ProjectInfoProjectName] = None
    project_description: str = None
    project_type: Optional[str] = None
    project_contributors: Optional[Union[str, list[str]]] = empty_list()
    project_collector_chief_scientist: Optional[str] = None
    BioProject_accession: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.project_name):
            self.MissingRequiredField("project_name")
        if not isinstance(self.project_name, ProjectInfoProjectName):
            self.project_name = ProjectInfoProjectName(self.project_name)

        if self._is_empty(self.project_description):
            self.MissingRequiredField("project_description")
        if not isinstance(self.project_description, str):
            self.project_description = str(self.project_description)

        if self.project_type is not None and not isinstance(self.project_type, str):
            self.project_type = str(self.project_type)

        if not isinstance(self.project_contributors, list):
            self.project_contributors = [self.project_contributors] if self.project_contributors is not None else []
        self.project_contributors = [v if isinstance(v, str) else str(v) for v in self.project_contributors]

        if self.project_collector_chief_scientist is not None and not isinstance(self.project_collector_chief_scientist, str):
            self.project_collector_chief_scientist = str(self.project_collector_chief_scientist)

        if self.BioProject_accession is not None and not isinstance(self.BioProject_accession, str):
            self.BioProject_accession = str(self.BioProject_accession)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class SpecimenInfo(YAMLRoot):
    """
    Information on specimen info
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["SpecimenInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:SpecimenInfo"
    class_name: ClassVar[str] = "SpecimenInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.SpecimenInfo

    specimen_name: Union[str, SpecimenInfoSpecimenName] = None
    specimen_taxon_id: Union[int, list[int]] = None
    host_taxon_id: int = None
    collection_date: str = None
    collection_country: str = None
    project_id: int = None
    host_subject_id: Optional[int] = None
    alternate_identifiers: Optional[Union[str, list[str]]] = empty_list()
    host_sex: Optional[str] = None
    specimen_accession: Optional[str] = None
    gravid: Optional[Union[bool, Bool]] = None
    blood_meal: Optional[Union[bool, Bool]] = None
    gravidity: Optional[int] = None
    microscopy_parasite_density_info: Optional[Union[Union[dict, ParasiteDensity], list[Union[dict, ParasiteDensity]]]] = empty_list()
    host_age: Optional[float] = None
    geo_admin1: Optional[str] = None
    geo_admin2: Optional[str] = None
    geo_admin3: Optional[str] = None
    lat_lon: Optional[str] = None
    specimen_store_loc: Optional[str] = None
    specimen_collect_device: Optional[str] = None
    specimen_type: Optional[str] = None
    specimen_comments: Optional[Union[str, list[str]]] = empty_list()
    travel_out_six_month: Optional[Union[str, list[str]]] = empty_list()
    storage_plate_info: Optional[Union[dict, PlateInfo]] = None
    env_medium: Optional[str] = None
    env_local_scale: Optional[str] = None
    env_broad_scale: Optional[str] = None
    drug_usage: Optional[Union[str, list[str]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.specimen_name):
            self.MissingRequiredField("specimen_name")
        if not isinstance(self.specimen_name, SpecimenInfoSpecimenName):
            self.specimen_name = SpecimenInfoSpecimenName(self.specimen_name)

        if self._is_empty(self.specimen_taxon_id):
            self.MissingRequiredField("specimen_taxon_id")
        if not isinstance(self.specimen_taxon_id, list):
            self.specimen_taxon_id = [self.specimen_taxon_id] if self.specimen_taxon_id is not None else []
        self.specimen_taxon_id = [v if isinstance(v, int) else int(v) for v in self.specimen_taxon_id]

        if self._is_empty(self.host_taxon_id):
            self.MissingRequiredField("host_taxon_id")
        if not isinstance(self.host_taxon_id, int):
            self.host_taxon_id = int(self.host_taxon_id)

        if self._is_empty(self.collection_date):
            self.MissingRequiredField("collection_date")
        if not isinstance(self.collection_date, str):
            self.collection_date = str(self.collection_date)

        if self._is_empty(self.collection_country):
            self.MissingRequiredField("collection_country")
        if not isinstance(self.collection_country, str):
            self.collection_country = str(self.collection_country)

        if self._is_empty(self.project_id):
            self.MissingRequiredField("project_id")
        if not isinstance(self.project_id, int):
            self.project_id = int(self.project_id)

        if self.host_subject_id is not None and not isinstance(self.host_subject_id, int):
            self.host_subject_id = int(self.host_subject_id)

        if not isinstance(self.alternate_identifiers, list):
            self.alternate_identifiers = [self.alternate_identifiers] if self.alternate_identifiers is not None else []
        self.alternate_identifiers = [v if isinstance(v, str) else str(v) for v in self.alternate_identifiers]

        if self.host_sex is not None and not isinstance(self.host_sex, str):
            self.host_sex = str(self.host_sex)

        if self.specimen_accession is not None and not isinstance(self.specimen_accession, str):
            self.specimen_accession = str(self.specimen_accession)

        if self.gravid is not None and not isinstance(self.gravid, Bool):
            self.gravid = Bool(self.gravid)

        if self.blood_meal is not None and not isinstance(self.blood_meal, Bool):
            self.blood_meal = Bool(self.blood_meal)

        if self.gravidity is not None and not isinstance(self.gravidity, int):
            self.gravidity = int(self.gravidity)

        if not isinstance(self.microscopy_parasite_density_info, list):
            self.microscopy_parasite_density_info = [self.microscopy_parasite_density_info] if self.microscopy_parasite_density_info is not None else []
        self.microscopy_parasite_density_info = [v if isinstance(v, ParasiteDensity) else ParasiteDensity(**as_dict(v)) for v in self.microscopy_parasite_density_info]

        if self.host_age is not None and not isinstance(self.host_age, float):
            self.host_age = float(self.host_age)

        if self.geo_admin1 is not None and not isinstance(self.geo_admin1, str):
            self.geo_admin1 = str(self.geo_admin1)

        if self.geo_admin2 is not None and not isinstance(self.geo_admin2, str):
            self.geo_admin2 = str(self.geo_admin2)

        if self.geo_admin3 is not None and not isinstance(self.geo_admin3, str):
            self.geo_admin3 = str(self.geo_admin3)

        if self.lat_lon is not None and not isinstance(self.lat_lon, str):
            self.lat_lon = str(self.lat_lon)

        if self.specimen_store_loc is not None and not isinstance(self.specimen_store_loc, str):
            self.specimen_store_loc = str(self.specimen_store_loc)

        if self.specimen_collect_device is not None and not isinstance(self.specimen_collect_device, str):
            self.specimen_collect_device = str(self.specimen_collect_device)

        if self.specimen_type is not None and not isinstance(self.specimen_type, str):
            self.specimen_type = str(self.specimen_type)

        if not isinstance(self.specimen_comments, list):
            self.specimen_comments = [self.specimen_comments] if self.specimen_comments is not None else []
        self.specimen_comments = [v if isinstance(v, str) else str(v) for v in self.specimen_comments]

        if not isinstance(self.travel_out_six_month, list):
            self.travel_out_six_month = [self.travel_out_six_month] if self.travel_out_six_month is not None else []
        self.travel_out_six_month = [v if isinstance(v, str) else str(v) for v in self.travel_out_six_month]

        if self.storage_plate_info is not None and not isinstance(self.storage_plate_info, PlateInfo):
            self.storage_plate_info = PlateInfo(**as_dict(self.storage_plate_info))

        if self.env_medium is not None and not isinstance(self.env_medium, str):
            self.env_medium = str(self.env_medium)

        if self.env_local_scale is not None and not isinstance(self.env_local_scale, str):
            self.env_local_scale = str(self.env_local_scale)

        if self.env_broad_scale is not None and not isinstance(self.env_broad_scale, str):
            self.env_broad_scale = str(self.env_broad_scale)

        if not isinstance(self.drug_usage, list):
            self.drug_usage = [self.drug_usage] if self.drug_usage is not None else []
        self.drug_usage = [v if isinstance(v, str) else str(v) for v in self.drug_usage]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class BioinformaticsRunInfo(YAMLRoot):
    """
    Information about the pipeline run that generated some of the microhaplotype detected and reads_by_stage
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["BioinformaticsRunInfo"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:BioinformaticsRunInfo"
    class_name: ClassVar[str] = "BioinformaticsRunInfo"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.BioinformaticsRunInfo

    bioinformatics_methods_id: int = None
    bioinformatics_run_name: str = None
    run_date: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.bioinformatics_methods_id):
            self.MissingRequiredField("bioinformatics_methods_id")
        if not isinstance(self.bioinformatics_methods_id, int):
            self.bioinformatics_methods_id = int(self.bioinformatics_methods_id)

        if self._is_empty(self.bioinformatics_run_name):
            self.MissingRequiredField("bioinformatics_run_name")
        if not isinstance(self.bioinformatics_run_name, str):
            self.bioinformatics_run_name = str(self.bioinformatics_run_name)

        if self.run_date is not None and not isinstance(self.run_date, str):
            self.run_date = str(self.run_date)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class PmoGenerationMethod(YAMLRoot):
    """
    Information about how a PMO was generated
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PmoGenerationMethod"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PmoGenerationMethod"
    class_name: ClassVar[str] = "PmoGenerationMethod"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PmoGenerationMethod

    program_version: str = None
    program_name: str = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.program_version):
            self.MissingRequiredField("program_version")
        if not isinstance(self.program_version, str):
            self.program_version = str(self.program_version)

        if self._is_empty(self.program_name):
            self.MissingRequiredField("program_name")
        if not isinstance(self.program_name, str):
            self.program_name = str(self.program_name)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class PmoHeader(YAMLRoot):
    """
    Information on the PMO file
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PmoHeader"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PmoHeader"
    class_name: ClassVar[str] = "PmoHeader"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PmoHeader

    pmo_version: str = None
    creation_date: Optional[str] = None
    generation_method: Optional[Union[dict, PmoGenerationMethod]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.pmo_version):
            self.MissingRequiredField("pmo_version")
        if not isinstance(self.pmo_version, str):
            self.pmo_version = str(self.pmo_version)

        if self.creation_date is not None and not isinstance(self.creation_date, str):
            self.creation_date = str(self.creation_date)

        if self.generation_method is not None and not isinstance(self.generation_method, PmoGenerationMethod):
            self.generation_method = PmoGenerationMethod(**as_dict(self.generation_method))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class StageReadCounts(YAMLRoot):
    """
    Information on the reads counts at several stages
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["StageReadCounts"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:StageReadCounts"
    class_name: ClassVar[str] = "StageReadCounts"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.StageReadCounts

    read_count: int = None
    stage: str = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.read_count):
            self.MissingRequiredField("read_count")
        if not isinstance(self.read_count, int):
            self.read_count = int(self.read_count)

        if self._is_empty(self.stage):
            self.MissingRequiredField("stage")
        if not isinstance(self.stage, str):
            self.stage = str(self.stage)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ReadCountsByStageForTarget(YAMLRoot):
    """
    Information on the reads counts at several stages of a pipeline for a target
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReadCountsByStageForTarget"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReadCountsByStageForTarget"
    class_name: ClassVar[str] = "ReadCountsByStageForTarget"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReadCountsByStageForTarget

    target_id: int = None
    stages: Union[Union[dict, StageReadCounts], list[Union[dict, StageReadCounts]]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
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


@dataclass(repr=False)
class ReadCountsByStageForLibrarySample(YAMLRoot):
    """
    Information on the reads counts at several stages of a pipeline for a library_sample
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReadCountsByStageForLibrarySample"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReadCountsByStageForLibrarySample"
    class_name: ClassVar[str] = "ReadCountsByStageForLibrarySample"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReadCountsByStageForLibrarySample

    library_sample_id: int = None
    total_raw_count: int = None
    read_counts_for_targets: Optional[Union[Union[dict, ReadCountsByStageForTarget], list[Union[dict, ReadCountsByStageForTarget]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.library_sample_id):
            self.MissingRequiredField("library_sample_id")
        if not isinstance(self.library_sample_id, int):
            self.library_sample_id = int(self.library_sample_id)

        if self._is_empty(self.total_raw_count):
            self.MissingRequiredField("total_raw_count")
        if not isinstance(self.total_raw_count, int):
            self.total_raw_count = int(self.total_raw_count)

        if not isinstance(self.read_counts_for_targets, list):
            self.read_counts_for_targets = [self.read_counts_for_targets] if self.read_counts_for_targets is not None else []
        self.read_counts_for_targets = [v if isinstance(v, ReadCountsByStageForTarget) else ReadCountsByStageForTarget(**as_dict(v)) for v in self.read_counts_for_targets]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ReadCountsByStage(YAMLRoot):
    """
    Information on the reads counts at several stages of a pipeline
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["ReadCountsByStage"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:ReadCountsByStage"
    class_name: ClassVar[str] = "ReadCountsByStage"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.ReadCountsByStage

    bioinformatics_run_id: int = None
    read_counts_by_library_sample_by_stage: Union[Union[dict, ReadCountsByStageForLibrarySample], list[Union[dict, ReadCountsByStageForLibrarySample]]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.bioinformatics_run_id):
            self.MissingRequiredField("bioinformatics_run_id")
        if not isinstance(self.bioinformatics_run_id, int):
            self.bioinformatics_run_id = int(self.bioinformatics_run_id)

        if self._is_empty(self.read_counts_by_library_sample_by_stage):
            self.MissingRequiredField("read_counts_by_library_sample_by_stage")
        if not isinstance(self.read_counts_by_library_sample_by_stage, list):
            self.read_counts_by_library_sample_by_stage = [self.read_counts_by_library_sample_by_stage] if self.read_counts_by_library_sample_by_stage is not None else []
        self.read_counts_by_library_sample_by_stage = [v if isinstance(v, ReadCountsByStageForLibrarySample) else ReadCountsByStageForLibrarySample(**as_dict(v)) for v in self.read_counts_by_library_sample_by_stage]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class PortableMicrohaplotypeObject(YAMLRoot):
    """
    Information on final results from a targeted amplicon analysis
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT["PortableMicrohaplotypeObject"]
    class_class_curie: ClassVar[str] = "portable_microhaplotype_object:PortableMicrohaplotypeObject"
    class_name: ClassVar[str] = "PortableMicrohaplotypeObject"
    class_model_uri: ClassVar[URIRef] = PORTABLE_MICROHAPLOTYPE_OBJECT.PortableMicrohaplotypeObject

    library_info: Union[dict[Union[str, LibrarySampleInfoLibrarySampleName], Union[dict, LibrarySampleInfo]], list[Union[dict, LibrarySampleInfo]]] = empty_dict()
    specimen_info: Union[dict[Union[str, SpecimenInfoSpecimenName], Union[dict, SpecimenInfo]], list[Union[dict, SpecimenInfo]]] = empty_dict()
    sequencing_info: Union[Union[dict, SequencingInfo], list[Union[dict, SequencingInfo]]] = None
    panel_info: Union[Union[dict, PanelInfo], list[Union[dict, PanelInfo]]] = None
    target_info: Union[Union[dict, TargetInfo], list[Union[dict, TargetInfo]]] = None
    targeted_genomes: Union[Union[dict, GenomeInfo], list[Union[dict, GenomeInfo]]] = None
    representative_microhaplotypes: Union[dict, RepresentativeMicrohaplotypes] = None
    bioinformatics_methods_info: Union[Union[dict, BioinformaticsMethodInfo], list[Union[dict, BioinformaticsMethodInfo]]] = None
    bioinformatics_run_info: Union[Union[dict, BioinformaticsRunInfo], list[Union[dict, BioinformaticsRunInfo]]] = None
    detected_microhaplotypes: Union[Union[dict, DetectedMicrohaplotypes], list[Union[dict, DetectedMicrohaplotypes]]] = None
    project_info: Union[dict[Union[str, ProjectInfoProjectName], Union[dict, ProjectInfo]], list[Union[dict, ProjectInfo]]] = empty_dict()
    pmo_header: Union[dict, PmoHeader] = None
    read_counts_by_stage: Optional[Union[Union[dict, ReadCountsByStage], list[Union[dict, ReadCountsByStage]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.library_info):
            self.MissingRequiredField("library_info")
        self._normalize_inlined_as_list(slot_name="library_info", slot_type=LibrarySampleInfo, key_name="library_sample_name", keyed=True)

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

        if self._is_empty(self.representative_microhaplotypes):
            self.MissingRequiredField("representative_microhaplotypes")
        if not isinstance(self.representative_microhaplotypes, RepresentativeMicrohaplotypes):
            self.representative_microhaplotypes = RepresentativeMicrohaplotypes(**as_dict(self.representative_microhaplotypes))

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

        if self._is_empty(self.detected_microhaplotypes):
            self.MissingRequiredField("detected_microhaplotypes")
        if not isinstance(self.detected_microhaplotypes, list):
            self.detected_microhaplotypes = [self.detected_microhaplotypes] if self.detected_microhaplotypes is not None else []
        self.detected_microhaplotypes = [v if isinstance(v, DetectedMicrohaplotypes) else DetectedMicrohaplotypes(**as_dict(v)) for v in self.detected_microhaplotypes]

        if self._is_empty(self.project_info):
            self.MissingRequiredField("project_info")
        self._normalize_inlined_as_list(slot_name="project_info", slot_type=ProjectInfo, key_name="project_name", keyed=True)

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

slots.library_sample_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_sample_id, name="library_sample_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_sample_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_sample_id, domain=None, range=int,
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

slots.seq = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq, name="seq", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq, domain=None, range=str,
                   pattern=re.compile(r'^[A-z]+$'))

slots.program_version = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_version, name="program_version", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program_version'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_version, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.markerOfInterest__marker_location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.marker_location, name="markerOfInterest__marker_location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('marker_location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.markerOfInterest__marker_location, domain=None, range=Union[dict, GenomicLocation])

slots.markerOfInterest__associations = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.associations, name="markerOfInterest__associations", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('associations'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.markerOfInterest__associations, domain=None, range=Optional[Union[str, list[str]]],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.targetInfo__target_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_name, name="targetInfo__target_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__target_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.targetInfo__gene_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gene_name, name="targetInfo__gene_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gene_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__gene_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.targetInfo__insert_location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.insert_location, name="targetInfo__insert_location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('insert_location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__insert_location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.targetInfo__forward_primer = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.forward_primer, name="targetInfo__forward_primer", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('forward_primer'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__forward_primer, domain=None, range=Union[dict, PrimerInfo])

slots.targetInfo__reverse_primer = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reverse_primer, name="targetInfo__reverse_primer", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reverse_primer'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__reverse_primer, domain=None, range=Union[dict, PrimerInfo])

slots.targetInfo__markers_of_interest = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.markers_of_interest, name="targetInfo__markers_of_interest", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('markers_of_interest'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__markers_of_interest, domain=None, range=Optional[Union[Union[dict, MarkerOfInterest], list[Union[dict, MarkerOfInterest]]]])

slots.targetInfo__target_attributes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_attributes, name="targetInfo__target_attributes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_attributes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targetInfo__target_attributes, domain=None, range=Optional[Union[str, list[str]]])

slots.reactionInfo__panel_targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_targets, name="reactionInfo__panel_targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reactionInfo__panel_targets, domain=None, range=Union[int, list[int]])

slots.reactionInfo__reaction_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reaction_name, name="reactionInfo__reaction_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reaction_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reactionInfo__reaction_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.panelInfo__reactions = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reactions, name="panelInfo__reactions", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reactions'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__reactions, domain=None, range=Union[Union[dict, ReactionInfo], list[Union[dict, ReactionInfo]]])

slots.panelInfo__panel_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_name, name="panelInfo__panel_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panelInfo__panel_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.pseudocigar__pseudocigar_seq = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar_seq, name="pseudocigar__pseudocigar_seq", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pseudocigar_seq'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar__pseudocigar_seq, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.pseudocigar__ref_loc = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.ref_loc, name="pseudocigar__ref_loc", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('ref_loc'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar__ref_loc, domain=None, range=Union[dict, GenomicLocation])

slots.pseudocigar__pseudocigar_generation_description = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar_generation_description, name="pseudocigar__pseudocigar_generation_description", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pseudocigar_generation_description'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar__pseudocigar_generation_description, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.representativeMicrohaplotype__microhaplotype_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotype_name, name="representativeMicrohaplotype__microhaplotype_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotype_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__microhaplotype_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.representativeMicrohaplotype__quality = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.quality, name="representativeMicrohaplotype__quality", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('quality'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__quality, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.representativeMicrohaplotype__pseudocigar = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pseudocigar, name="representativeMicrohaplotype__pseudocigar", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pseudocigar'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__pseudocigar, domain=None, range=Optional[Union[dict, Pseudocigar]],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.representativeMicrohaplotype__masking = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.masking, name="representativeMicrohaplotype__masking", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('masking'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__masking, domain=None, range=Optional[Union[Union[dict, MaskingInfo], list[Union[dict, MaskingInfo]]]])

slots.representativeMicrohaplotype__alt_annotations = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alt_annotations, name="representativeMicrohaplotype__alt_annotations", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alt_annotations'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotype__alt_annotations, domain=None, range=Optional[Union[str, list[str]]])

slots.maskingInfo__seq_start = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_start, name="maskingInfo__seq_start", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_start'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__seq_start, domain=None, range=int)

slots.maskingInfo__seq_segment_size = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_segment_size, name="maskingInfo__seq_segment_size", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_segment_size'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__seq_segment_size, domain=None, range=int)

slots.maskingInfo__replacement_size = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.replacement_size, name="maskingInfo__replacement_size", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('replacement_size'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__replacement_size, domain=None, range=int)

slots.maskingInfo__masking_generation_description = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.masking_generation_description, name="maskingInfo__masking_generation_description", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('masking_generation_description'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.maskingInfo__masking_generation_description, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.representativeMicrohaplotypes__targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targets, name="representativeMicrohaplotypes__targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypes__targets, domain=None, range=Union[dict[Union[int, RepresentativeMicrohaplotypesForTargetTargetId], Union[dict, RepresentativeMicrohaplotypesForTarget]], list[Union[dict, RepresentativeMicrohaplotypesForTarget]]])

slots.representativeMicrohaplotypesForTarget__microhaplotypes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypes, name="representativeMicrohaplotypesForTarget__microhaplotypes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microhaplotypes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypesForTarget__microhaplotypes, domain=None, range=Union[Union[dict, RepresentativeMicrohaplotype], list[Union[dict, RepresentativeMicrohaplotype]]])

slots.representativeMicrohaplotypesForTarget__mhap_location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.mhap_location, name="representativeMicrohaplotypesForTarget__mhap_location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('mhap_location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representativeMicrohaplotypesForTarget__mhap_location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.detectedMicrohaplotypes__library_samples = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_samples, name="detectedMicrohaplotypes__library_samples", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_samples'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.detectedMicrohaplotypes__library_samples, domain=None, range=Union[Union[dict, DetectedMicrohaplotypesForSample], list[Union[dict, DetectedMicrohaplotypesForSample]]])

slots.genomeInfo__name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.name, name="genomeInfo__name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.genomeInfo__genome_version = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genome_version, name="genomeInfo__genome_version", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('genome_version'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__genome_version, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.genomeInfo__taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.taxon_id, name="genomeInfo__taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__taxon_id, domain=None, range=Union[int, list[int]],
                   pattern=re.compile(r'^[0-9]$'))

slots.genomeInfo__url = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.url, name="genomeInfo__url", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('url'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__url, domain=None, range=str,
                   pattern=re.compile(r'^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$'))

slots.genomeInfo__chromosomes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.chromosomes, name="genomeInfo__chromosomes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('chromosomes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomeInfo__chromosomes, domain=None, range=Optional[Union[str, list[str]]],
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

slots.genomicLocation__alt_seq = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alt_seq, name="genomicLocation__alt_seq", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alt_seq'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.genomicLocation__alt_seq, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-]+$'))

slots.primerInfo__location = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.location, name="primerInfo__location", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('location'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.primerInfo__location, domain=None, range=Optional[Union[dict, GenomicLocation]])

slots.detectedMicrohaplotypesForSample__target_results = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_results, name="detectedMicrohaplotypesForSample__target_results", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_results'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.detectedMicrohaplotypesForSample__target_results, domain=None, range=Union[Union[dict, DetectedMicrohaplotypesForTarget], list[Union[dict, DetectedMicrohaplotypesForTarget]]])

slots.microhaplotypeForTarget__reads = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.reads, name="microhaplotypeForTarget__reads", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('reads'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__reads, domain=None, range=int)

slots.microhaplotypeForTarget__umis = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.umis, name="microhaplotypeForTarget__umis", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('umis'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microhaplotypeForTarget__umis, domain=None, range=Optional[int])

slots.detectedMicrohaplotypesForTarget__mhaps = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.mhaps, name="detectedMicrohaplotypesForTarget__mhaps", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('mhaps'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.detectedMicrohaplotypesForTarget__mhaps, domain=None, range=Union[Union[dict, MicrohaplotypeForTarget], list[Union[dict, MicrohaplotypeForTarget]]])

slots.bioinformaticsMethodInfo__demultiplexing_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.demultiplexing_method, name="bioinformaticsMethodInfo__demultiplexing_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('demultiplexing_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__demultiplexing_method, domain=None, range=Union[dict, BioMethod])

slots.bioinformaticsMethodInfo__denoising_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.denoising_method, name="bioinformaticsMethodInfo__denoising_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('denoising_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__denoising_method, domain=None, range=Union[dict, BioMethod])

slots.bioinformaticsMethodInfo__additional_methods = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.additional_methods, name="bioinformaticsMethodInfo__additional_methods", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('additional_methods'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__additional_methods, domain=None, range=Optional[Union[Union[dict, BioMethod], list[Union[dict, BioMethod]]]])

slots.bioinformaticsMethodInfo__bioinformatics_method_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_method_name, name="bioinformaticsMethodInfo__bioinformatics_method_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_method_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsMethodInfo__bioinformatics_method_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.bioMethod__program = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program, name="bioMethod__program", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__program, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.bioMethod__program_description = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_description, name="bioMethod__program_description", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program_description'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__program_description, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.bioMethod__program_url = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.program_url, name="bioMethod__program_url", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('program_url'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__program_url, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.bioMethod__additional_argument = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.additional_argument, name="bioMethod__additional_argument", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('additional_argument'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioMethod__additional_argument, domain=None, range=Optional[Union[str, list[str]]],
                   pattern=re.compile(r'^[A-z-._0-9{}\(\),\/\ ]+$'))

slots.plateInfo__plate_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_name, name="plateInfo__plate_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plateInfo__plate_name, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.plateInfo__plate_row = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_row, name="plateInfo__plate_row", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_row'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plateInfo__plate_row, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z]$'))

slots.plateInfo__plate_col = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plate_col, name="plateInfo__plate_col", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('plate_col'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.plateInfo__plate_col, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]+$'))

slots.librarySampleInfo__fastqs_loc = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.fastqs_loc, name="librarySampleInfo__fastqs_loc", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('fastqs_loc'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.librarySampleInfo__fastqs_loc, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-.;,_0-9\(\),\/\ ]+$'))

slots.librarySampleInfo__run_accession = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.run_accession, name="librarySampleInfo__run_accession", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('run_accession'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.librarySampleInfo__run_accession, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9]+$'))

slots.librarySampleInfo__library_sample_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_sample_name, name="librarySampleInfo__library_sample_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_sample_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.librarySampleInfo__library_sample_name, domain=None, range=URIRef,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.librarySampleInfo__library_prep_plate_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_prep_plate_info, name="librarySampleInfo__library_prep_plate_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_prep_plate_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.librarySampleInfo__library_prep_plate_info, domain=None, range=Optional[Union[dict, PlateInfo]])

slots.librarySampleInfo__qpcr_parasite_density_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.qpcr_parasite_density_info, name="librarySampleInfo__qpcr_parasite_density_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('qpcr_parasite_density_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.librarySampleInfo__qpcr_parasite_density_info, domain=None, range=Optional[Union[Union[dict, ParasiteDensity], list[Union[dict, ParasiteDensity]]]])

slots.sequencingInfo__sequencing_info_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info_name, name="sequencingInfo__sequencing_info_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__sequencing_info_name, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__seq_platform = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_platform, name="sequencingInfo__seq_platform", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_platform'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_platform, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__seq_instrument_model = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_instrument_model, name="sequencingInfo__seq_instrument_model", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_instrument_model'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_instrument_model, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__seq_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_date, name="sequencingInfo__seq_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_date, domain=None, range=Optional[str],
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

slots.sequencingInfo__library_screen = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_screen, name="sequencingInfo__library_screen", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_screen'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__library_screen, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.sequencingInfo__library_kit = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_kit, name="sequencingInfo__library_kit", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_kit'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__library_kit, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.sequencingInfo__library_layout = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_layout, name="sequencingInfo__library_layout", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_layout'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__library_layout, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__library_strategy = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_strategy, name="sequencingInfo__library_strategy", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_strategy'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__library_strategy, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__library_source = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_source, name="sequencingInfo__library_source", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_source'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__library_source, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__library_selection = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_selection, name="sequencingInfo__library_selection", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_selection'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__library_selection, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.sequencingInfo__seq_center = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.seq_center, name="sequencingInfo__seq_center", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('seq_center'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencingInfo__seq_center, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.parasiteDensity__density_method = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.density_method, name="parasiteDensity__density_method", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('density_method'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasiteDensity__density_method, domain=None, range=str,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.parasiteDensity__parasite_density = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasite_density, name="parasiteDensity__parasite_density", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('parasite_density'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasiteDensity__parasite_density, domain=None, range=float,
                   pattern=re.compile(r'^[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?$'))

slots.parasiteDensity__date_measured = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.date_measured, name="parasiteDensity__date_measured", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('date_measured'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasiteDensity__date_measured, domain=None, range=Optional[str],
                   pattern=re.compile(r'(?:\d{4}(?:-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?)?|NA)'))

slots.parasiteDensity__density_method_comments = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.density_method_comments, name="parasiteDensity__density_method_comments", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('density_method_comments'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.parasiteDensity__density_method_comments, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.projectInfo__project_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_name, name="projectInfo__project_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.projectInfo__project_name, domain=None, range=URIRef,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.projectInfo__project_description = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_description, name="projectInfo__project_description", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_description'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.projectInfo__project_description, domain=None, range=str)

slots.projectInfo__project_type = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_type, name="projectInfo__project_type", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_type'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.projectInfo__project_type, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.projectInfo__project_contributors = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_contributors, name="projectInfo__project_contributors", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_contributors'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.projectInfo__project_contributors, domain=None, range=Optional[Union[str, list[str]]],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.projectInfo__project_collector_chief_scientist = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_collector_chief_scientist, name="projectInfo__project_collector_chief_scientist", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_collector_chief_scientist'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.projectInfo__project_collector_chief_scientist, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9;|\(\),\/\ ]+$'))

slots.projectInfo__BioProject_accession = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.BioProject_accession, name="projectInfo__BioProject_accession", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('BioProject_accession'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.projectInfo__BioProject_accession, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__specimen_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_name, name="specimenInfo__specimen_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_name, domain=None, range=URIRef,
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__specimen_taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_taxon_id, name="specimenInfo__specimen_taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_taxon_id, domain=None, range=Union[int, list[int]],
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__host_subject_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.host_subject_id, name="specimenInfo__host_subject_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('host_subject_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__host_subject_id, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__host_taxon_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.host_taxon_id, name="specimenInfo__host_taxon_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('host_taxon_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__host_taxon_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__alternate_identifiers = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.alternate_identifiers, name="specimenInfo__alternate_identifiers", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('alternate_identifiers'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__alternate_identifiers, domain=None, range=Optional[Union[str, list[str]]],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__host_sex = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.host_sex, name="specimenInfo__host_sex", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('host_sex'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__host_sex, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__specimen_accession = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_accession, name="specimenInfo__specimen_accession", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_accession'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_accession, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9 ]+$'))

slots.specimenInfo__gravid = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gravid, name="specimenInfo__gravid", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gravid'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__gravid, domain=None, range=Optional[Union[bool, Bool]])

slots.specimenInfo__blood_meal = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.blood_meal, name="specimenInfo__blood_meal", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('blood_meal'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__blood_meal, domain=None, range=Optional[Union[bool, Bool]])

slots.specimenInfo__gravidity = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.gravidity, name="specimenInfo__gravidity", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('gravidity'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__gravidity, domain=None, range=Optional[int],
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__microscopy_parasite_density_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.microscopy_parasite_density_info, name="specimenInfo__microscopy_parasite_density_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('microscopy_parasite_density_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__microscopy_parasite_density_info, domain=None, range=Optional[Union[Union[dict, ParasiteDensity], list[Union[dict, ParasiteDensity]]]])

slots.specimenInfo__collection_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collection_date, name="specimenInfo__collection_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collection_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collection_date, domain=None, range=str,
                   pattern=re.compile(r'(?:\d{4}(?:-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?)?|NA)'))

slots.specimenInfo__host_age = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.host_age, name="specimenInfo__host_age", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('host_age'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__host_age, domain=None, range=Optional[float],
                   pattern=re.compile(r'^\d*\.?\d+$'))

slots.specimenInfo__collection_country = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.collection_country, name="specimenInfo__collection_country", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('collection_country'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__collection_country, domain=None, range=str,
                   pattern=re.compile(r'^[A-Za-z0-9 ,._:'–-]+$'))

slots.specimenInfo__geo_admin1 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin1, name="specimenInfo__geo_admin1", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin1'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin1, domain=None, range=Optional[str])

slots.specimenInfo__geo_admin2 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin2, name="specimenInfo__geo_admin2", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin2'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin2, domain=None, range=Optional[str])

slots.specimenInfo__geo_admin3 = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.geo_admin3, name="specimenInfo__geo_admin3", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('geo_admin3'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__geo_admin3, domain=None, range=Optional[str])

slots.specimenInfo__lat_lon = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.lat_lon, name="specimenInfo__lat_lon", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('lat_lon'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__lat_lon, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[-+]?\d{1,2}(?:\.\d+)?,[-+]?\d{1,3}(?:\.\d+)?$'))

slots.specimenInfo__specimen_store_loc = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_store_loc, name="specimenInfo__specimen_store_loc", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_store_loc'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_store_loc, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__specimen_collect_device = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_collect_device, name="specimenInfo__specimen_collect_device", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_collect_device'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_collect_device, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__specimen_type = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_type, name="specimenInfo__specimen_type", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_type'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_type, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__project_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_id, name="specimenInfo__project_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__project_id, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.specimenInfo__specimen_comments = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_comments, name="specimenInfo__specimen_comments", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_comments'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__specimen_comments, domain=None, range=Optional[Union[str, list[str]]],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__travel_out_six_month = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.travel_out_six_month, name="specimenInfo__travel_out_six_month", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('travel_out_six_month'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__travel_out_six_month, domain=None, range=Optional[Union[str, list[str]]],
                   pattern=re.compile(r'^[A-z-._0-9\(\),\/\ ]+$'))

slots.specimenInfo__storage_plate_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.storage_plate_info, name="specimenInfo__storage_plate_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('storage_plate_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__storage_plate_info, domain=None, range=Optional[Union[dict, PlateInfo]])

slots.specimenInfo__env_medium = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.env_medium, name="specimenInfo__env_medium", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('env_medium'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__env_medium, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9;|\(\),\/\ ]+$'))

slots.specimenInfo__env_local_scale = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.env_local_scale, name="specimenInfo__env_local_scale", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('env_local_scale'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__env_local_scale, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9;|\(\),\/\ ]+$'))

slots.specimenInfo__env_broad_scale = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.env_broad_scale, name="specimenInfo__env_broad_scale", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('env_broad_scale'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__env_broad_scale, domain=None, range=Optional[str],
                   pattern=re.compile(r'^[A-z-._0-9;|\(\),\/\ ]+$'))

slots.specimenInfo__drug_usage = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.drug_usage, name="specimenInfo__drug_usage", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('drug_usage'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimenInfo__drug_usage, domain=None, range=Optional[Union[str, list[str]]],
                   pattern=re.compile(r'^[A-z-._0-9;|\(\),\/\ ]+$'))

slots.bioinformaticsRunInfo__run_date = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.run_date, name="bioinformaticsRunInfo__run_date", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('run_date'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsRunInfo__run_date, domain=None, range=Optional[str],
                   pattern=re.compile(r'\d{4}-(?:0[1-9]|1[0-2])(?:-(?:0[1-9]|[12][0-9]|3[01]))?'))

slots.bioinformaticsRunInfo__bioinformatics_run_name = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_run_name, name="bioinformaticsRunInfo__bioinformatics_run_name", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_run_name'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformaticsRunInfo__bioinformatics_run_name, domain=None, range=str,
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
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStageForTarget__stages, domain=None, range=Union[Union[dict, StageReadCounts], list[Union[dict, StageReadCounts]]])

slots.readCountsByStageForLibrarySample__total_raw_count = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.total_raw_count, name="readCountsByStageForLibrarySample__total_raw_count", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('total_raw_count'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStageForLibrarySample__total_raw_count, domain=None, range=int,
                   pattern=re.compile(r'^[0-9]+$'))

slots.readCountsByStageForLibrarySample__read_counts_for_targets = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_counts_for_targets, name="readCountsByStageForLibrarySample__read_counts_for_targets", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_counts_for_targets'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStageForLibrarySample__read_counts_for_targets, domain=None, range=Optional[Union[Union[dict, ReadCountsByStageForTarget], list[Union[dict, ReadCountsByStageForTarget]]]])

slots.readCountsByStage__read_counts_by_library_sample_by_stage = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_counts_by_library_sample_by_stage, name="readCountsByStage__read_counts_by_library_sample_by_stage", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_counts_by_library_sample_by_stage'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.readCountsByStage__read_counts_by_library_sample_by_stage, domain=None, range=Union[Union[dict, ReadCountsByStageForLibrarySample], list[Union[dict, ReadCountsByStageForLibrarySample]]])

slots.portableMicrohaplotypeObject__library_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.library_info, name="portableMicrohaplotypeObject__library_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('library_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__library_info, domain=None, range=Union[dict[Union[str, LibrarySampleInfoLibrarySampleName], Union[dict, LibrarySampleInfo]], list[Union[dict, LibrarySampleInfo]]])

slots.portableMicrohaplotypeObject__specimen_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_info, name="portableMicrohaplotypeObject__specimen_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__specimen_info, domain=None, range=Union[dict[Union[str, SpecimenInfoSpecimenName], Union[dict, SpecimenInfo]], list[Union[dict, SpecimenInfo]]])

slots.portableMicrohaplotypeObject__sequencing_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.sequencing_info, name="portableMicrohaplotypeObject__sequencing_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('sequencing_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__sequencing_info, domain=None, range=Union[Union[dict, SequencingInfo], list[Union[dict, SequencingInfo]]])

slots.portableMicrohaplotypeObject__panel_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.panel_info, name="portableMicrohaplotypeObject__panel_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('panel_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__panel_info, domain=None, range=Union[Union[dict, PanelInfo], list[Union[dict, PanelInfo]]])

slots.portableMicrohaplotypeObject__target_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_info, name="portableMicrohaplotypeObject__target_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__target_info, domain=None, range=Union[Union[dict, TargetInfo], list[Union[dict, TargetInfo]]])

slots.portableMicrohaplotypeObject__targeted_genomes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.targeted_genomes, name="portableMicrohaplotypeObject__targeted_genomes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('targeted_genomes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__targeted_genomes, domain=None, range=Union[Union[dict, GenomeInfo], list[Union[dict, GenomeInfo]]])

slots.portableMicrohaplotypeObject__representative_microhaplotypes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.representative_microhaplotypes, name="portableMicrohaplotypeObject__representative_microhaplotypes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('representative_microhaplotypes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__representative_microhaplotypes, domain=None, range=Union[dict, RepresentativeMicrohaplotypes])

slots.portableMicrohaplotypeObject__bioinformatics_methods_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_methods_info, name="portableMicrohaplotypeObject__bioinformatics_methods_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_methods_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__bioinformatics_methods_info, domain=None, range=Union[Union[dict, BioinformaticsMethodInfo], list[Union[dict, BioinformaticsMethodInfo]]])

slots.portableMicrohaplotypeObject__bioinformatics_run_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.bioinformatics_run_info, name="portableMicrohaplotypeObject__bioinformatics_run_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('bioinformatics_run_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__bioinformatics_run_info, domain=None, range=Union[Union[dict, BioinformaticsRunInfo], list[Union[dict, BioinformaticsRunInfo]]])

slots.portableMicrohaplotypeObject__detected_microhaplotypes = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.detected_microhaplotypes, name="portableMicrohaplotypeObject__detected_microhaplotypes", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('detected_microhaplotypes'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__detected_microhaplotypes, domain=None, range=Union[Union[dict, DetectedMicrohaplotypes], list[Union[dict, DetectedMicrohaplotypes]]])

slots.portableMicrohaplotypeObject__project_info = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.project_info, name="portableMicrohaplotypeObject__project_info", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('project_info'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__project_info, domain=None, range=Union[dict[Union[str, ProjectInfoProjectName], Union[dict, ProjectInfo]], list[Union[dict, ProjectInfo]]])

slots.portableMicrohaplotypeObject__pmo_header = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.pmo_header, name="portableMicrohaplotypeObject__pmo_header", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('pmo_header'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__pmo_header, domain=None, range=Union[dict, PmoHeader])

slots.portableMicrohaplotypeObject__read_counts_by_stage = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.read_counts_by_stage, name="portableMicrohaplotypeObject__read_counts_by_stage", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('read_counts_by_stage'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.portableMicrohaplotypeObject__read_counts_by_stage, domain=None, range=Optional[Union[Union[dict, ReadCountsByStage], list[Union[dict, ReadCountsByStage]]]])

slots.RepresentativeMicrohaplotypesForTarget_target_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.target_id, name="RepresentativeMicrohaplotypesForTarget_target_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('target_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.RepresentativeMicrohaplotypesForTarget_target_id, domain=RepresentativeMicrohaplotypesForTarget, range=Union[int, RepresentativeMicrohaplotypesForTargetTargetId],
                   pattern=re.compile(r'^[0-9]+$'))

slots.LibrarySampleInfo_specimen_id = Slot(uri=PORTABLE_MICROHAPLOTYPE_OBJECT.specimen_id, name="LibrarySampleInfo_specimen_id", curie=PORTABLE_MICROHAPLOTYPE_OBJECT.curie('specimen_id'),
                   model_uri=PORTABLE_MICROHAPLOTYPE_OBJECT.LibrarySampleInfo_specimen_id, domain=LibrarySampleInfo, range=int,
                   pattern=re.compile(r'^[0-9]+$'))