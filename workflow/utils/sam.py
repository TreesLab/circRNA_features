import logging
import re
from collections import defaultdict, namedtuple


logger = logging.getLogger(__name__)


class CIGAR:
    """Use to parse the CIGAR string."""

    _CIGAR_PAT = re.compile(r'(([0-9]+)([MIDNSHP=X]))')

    def __init__(self, cigar_str):
        self._cigar_str = cigar_str
        self._parse()

    def _parse(self):
        logger.debug(self._cigar_str)

        if self._cigar_str == '*':
            self._cigar_list = []
            self._cigar_dict = {}
        else:
            parsed_result = re.findall(self._CIGAR_PAT, self._cigar_str)
            logger.debug(parsed_result)

            self._cigar_list = [cigar for cigar, _, _ in parsed_result]

            self._cigar_dict = defaultdict(list)
            for _, num, op in parsed_result:
                self._cigar_dict[op].append(int(num))

        logger.debug(self._cigar_list)
        logger.debug(self._cigar_dict)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._cigar_list[key]

        elif isinstance(key, str):
            return self._cigar_dict.get(key, [])

        else:
            pass

    def __repr__(self):
        return str(self._cigar_list)

    def __str__(self):
        return self._cigar_str


class SamFormat:
    _optional_field_pattern = re.compile(r'([A-Za-z][A-Za-z0-9]):([AifZHB]):(.+)')
    _optional_field_tuple = namedtuple('OptField', ['tag', 'type_', 'value'])

    def __init__(self, sam_string):
        self._sam_string = sam_string.rstrip('\n')
        self._init()
        self._parse()

    def _init(self):
        self.is_header = False

        self.qname = None
        self.flag = None
        self.rname = None
        self.pos = None
        self.mapq = None
        self.cigar = None
        self.rnext = None
        self.pnext = None
        self.tlen = None
        self.seq = None
        self.qual = None
        self.optional_fields = None

    def _parse(self):
        if self._sam_string.startswith('@'):
            self.is_header = True
        else:
            data = self._sam_string.split('\t')

            self.qname = data[0]
            self.flag = int(data[1])
            self.rname = data[2]
            self.pos = int(data[3])
            self.mapq = int(data[4])
            self.cigar = CIGAR(data[5])
            self.rnext = data[6]
            self.pnext = int(data[7])
            self.tlen = int(data[8])
            self.seq = data[9]
            self.qual = data[10]

            self.optional_fields = self._parse_optional_fields(data[11:])

    def _parse_optional_fields(self, fields):
        fields_dict = {}
        for field in fields:
            m = re.search(self._optional_field_pattern, field)
            if m:
                tag, type_, value = m.groups()
                fields_dict[tag] = self._optional_field_tuple(tag, type_, value)

        return fields_dict

    def __repr__(self):
        return self._sam_string

    def __str__(self):
        return self._sam_string

    @property
    def is_unmapped(self):
        return str(self.cigar) == '*'
