import gzip
import hashlib
import re
from abc import ABC
from typing import List, Optional

from janis_core import File
from janis_unix import Gunzipped

from janis_bioinformatics.data_types.tabix import FileTabix
from janis_core.tool.test_classes import (
    TTestPreprocessor,
    TTestExpectedOutput,
    TTestCase,
)
import operator


class Vcf(File):
    def __init__(self, optional=False):
        super().__init__(optional=optional, extension=".vcf")

    @staticmethod
    def name():
        return "VCF"

    def doc(self):
        return """
    Variant Call Format:

    The Variant Call Format (VCF) specifies the format of a text file 
    used in bioinformatics for storing gene sequence variations. 

    Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf
    """.strip()

    # @classmethod
    # def md5(cls, file_path: str):
    #     with open(file_path, "r") as f:
    #         meaningful_content = False
    #         hash_md5 = hashlib.md5()
    #         while True:
    #             line = f.readline()
    #             if not line:
    #                 break
    #             if re.match("^#CHROM", line):
    #                 meaningful_content = True
    #             if not meaningful_content:
    #                 continue
    #             hash_md5.update(line.encode())
    #     return hash_md5.hexdigest()

    @classmethod
    def md5_without_header(cls, file_path: str, headers_to_remove: List[str]) -> str:
        """
        Compute md5 of a vcf file with unwanted headers removed

        :param file_path: path to the file
        :type file_path: str
        :param headers_to_remove: headers to be removed before computing md5
        :type headers_to_remove: List[str]
        :return: md5
        :rtype: str
        """
        with open(file_path, "r") as f:
            hash_md5 = hashlib.md5()
            while True:
                line = f.readline()
                if not line:
                    break
                if all(("##" + header) not in line for header in headers_to_remove):
                    hash_md5.update(line.encode())
        return hash_md5.hexdigest()

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_size: int,
        line_count: int,
        headers_to_remove: Optional[List[str]] = None,
        md5_value: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_size,
            ),
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.LineCount,
                operator=operator.eq,
                expected_value=line_count,
            ),
        ]
        if (headers_to_remove is not None) and (md5_value is not None):
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=Vcf.md5_without_header,
                    operator=operator.eq,
                    expected_value=md5_value,
                    additional_paras={"headers_to_remove": headers_to_remove},
                ),
            ]
        return outcome


class VcfIdx(Vcf):
    @staticmethod
    def name():
        return "IndexedVCF"

    @staticmethod
    def secondary_files():
        return [".idx"]

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_vcf_size: int,
        line_count: int,
        min_idx_size: int,
        headers_to_remove: Optional[List[str]] = None,
        vcf_md5: Optional[str] = None,
        idx_md5: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = super().basic_test(
            tag, min_vcf_size, line_count, headers_to_remove, vcf_md5
        ) + [
            TTestExpectedOutput(
                tag=tag,
                suffix_secondary_file=".idx",
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_idx_size,
            ),
        ]
        if idx_md5 is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileMd5,
                    operator=operator.eq,
                    expected_value=idx_md5,
                ),
            ]
        return outcome


class CompressedVcf(Gunzipped):
    def __init__(self, optional=False):
        super().__init__(inner_type=Vcf, optional=optional, extension=".vcf.gz")

    @staticmethod
    def name():
        return "CompressedVCF"

    def doc(self):
        return ".vcf.gz"

    # @classmethod
    # def md5(cls, file_path: str):
    #     with gzip.open(file_path, "rt") as f:
    #         meaningful_content = False
    #         hash_md5 = hashlib.md5()
    #         while True:
    #             line = f.readline()
    #             if not line:
    #                 break
    #             if re.match("^#CHROM", line):
    #                 meaningful_content = True
    #             if not meaningful_content:
    #                 continue
    #             hash_md5.update(line.encode())
    #     return hash_md5.hexdigest()

    @classmethod
    def LineCount(cls, file_path: str):
        count = 0
        with gzip.open(file_path, "rt") as f:
            while f.readline():
                count += 1
        return count

    @classmethod
    def md5_without_header(cls, file_path: str, headers_to_remove: List[str]) -> str:
        """
        Compute md5 of a vcf.gz file with unwanted headers removed

        :param file_path: path to the file
        :type file_path: str
        :param headers_to_remove: headers to be removed before computing md5
        :type headers_to_remove: List[str]
        :return: md5
        :rtype: str
        """
        with gzip.open(file_path, "rt") as f:
            hash_md5 = hashlib.md5()
            while True:
                line = f.readline()
                if not line:
                    break
                if all(("##" + header) not in line for header in headers_to_remove):
                    hash_md5.update(line.encode())
        return hash_md5.hexdigest()

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_size: int,
        line_count: int,
        headers_to_remove: Optional[List[str]] = None,
        md5_value: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_size,
            ),
            TTestExpectedOutput(
                tag=tag,
                preprocessor=CompressedVcf.LineCount,
                operator=operator.eq,
                expected_value=line_count,
            ),
        ]
        if (headers_to_remove is not None) and (md5_value is not None):
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=CompressedVcf.md5_without_header,
                    operator=operator.eq,
                    expected_value=md5_value,
                    additional_paras={"headers_to_remove": headers_to_remove},
                ),
            ]
        return outcome


class VcfTabix(CompressedVcf, FileTabix):
    @staticmethod
    def name():
        return "CompressedIndexedVCF"

    def doc(self):
        return ".vcf.gz with .vcf.gz.tbi file"

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_vcf_size: int,
        line_count: int,
        min_tbi_size: int,
        headers_to_remove: Optional[List[str]] = None,
        vcf_md5: Optional[str] = None,
        tbi_md5: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = super().basic_test(
            tag, min_vcf_size, line_count, headers_to_remove, vcf_md5
        ) + [
            TTestExpectedOutput(
                tag=tag,
                suffix_secondary_file=".tbi",
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_tbi_size,
            ),
        ]
        if tbi_md5 is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileMd5,
                    operator=operator.eq,
                    expected_value=tbi_md5,
                ),
            ]
        return outcome


# class GVCF(Vcf):
#     @staticmethod
#     def name():
#         return "gVCF"
#
#     def doc(self):
#         return """
# Section 5.5: Representing unspecified alleles and REF only blocks (gVCF)
# Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf
#         """
