#include "align_utility.h"

namespace is {

std::string CigarToString(const std::vector<CigarOp> &cigar) {
	std::stringstream ss;
	for (uint32_t i=0; i<cigar.size(); i++) {
		ss << cigar[i].count << ALN_OP_TO_CHAR[cigar[i].op]	;
	}
	return ss.str();
}

std::string CigarToBasicString(const std::vector<CigarOp> &cigar) {
	std::stringstream ss;
	std::vector<CigarOp> cigar_basic(cigar);

	int32_t offset = 0;
	for (int32_t i=1; i<(int32_t) cigar_basic.size(); i++) {
		if (ALN_OP_TO_BASIC_CHAR[cigar_basic[i].op] == ALN_OP_TO_BASIC_CHAR[cigar_basic[i-offset-1].op]) {
			cigar_basic[i-offset-1].count += cigar_basic[i].count;
			offset += 1;
			continue;
		}
		cigar_basic[i-offset] = cigar_basic[i];
	}
	for (int32_t i=0; i<(((int32_t) cigar_basic.size())-offset); i++) {
		ss << cigar_basic[i].count << ALN_OP_TO_BASIC_CHAR[cigar_basic[i].op]	;
	}
	return ss.str();
}

void CigarToEdlibAln(const std::vector<CigarOp> &cigar, std::vector<int8_t>& alignment) {
	alignment.clear();
	int32_t aln_len = 0;
	for (uint32_t i=0; i<cigar.size(); i++) {
		aln_len += cigar[i].count;
	}
	alignment.reserve(aln_len);
	for (uint32_t i=0; i<cigar.size(); i++) {
		std::vector<int8_t> temp(cigar[i].count, ALN_OP_TO_EDLIB[cigar[i].op]);
		alignment.insert(alignment.end(), temp.begin(), temp.end());
	}
}

void CigarToAlignment(const char* q, int64_t ql, int64_t q_start, int64_t q_end,
					  const char* t, int64_t tl, int64_t t_start, int64_t t_end, AlignType aln_type,
					  const std::vector<CigarOp> &cigar, std::string &alnq, std::string &alnt, std::string &alnm) {
	std::stringstream ssq, sst, ssm;

	int32_t posq = 0, post = 0;

  if (aln_type == kGlobal) {
    for (int32_t i=0; i<t_start; i++) {
      ssq << "-"; sst << t[post++]; ssm << " ";
    }
    for (int32_t i=0; i<q_start; i++) {
      ssq << q[posq++]; sst << "-"; ssm << " ";
    }
  }
  posq = q_start;
  post = t_start;

	for (uint32_t i=0; i<cigar.size(); i++) {
		if (cigar[i].op == ALN_OP_EQ || cigar[i].op == ALN_OP_X) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << q[posq++];	sst << t[post++];	ssm << ALN_OP_TO_MATCH[cigar[i].op];
			}
		}
		if (cigar[i].op == ALN_OP_I) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << q[posq++];	sst << "-";	ssm << " ";
			}
		}
		if (cigar[i].op == ALN_OP_D) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << "-";	sst << t[post++];	ssm << " ";
			}
		}
		if (cigar[i].op == ALN_OP_S) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << q[posq++];	sst << "-";	ssm << "-";
			}
		}

	}

  if (aln_type == kGlobal) {
    for (; post < tl; post++) {
      ssq << "-"; sst << t[post]; ssm << " ";
    }
    for (; posq < ql; posq++) {
      ssq << q[posq]; sst << "-"; ssm << " ";
    }
  }

	alnq = ssq.str();
	alnt = sst.str();
	alnm = ssm.str();
}

}
