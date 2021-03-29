#include <boost/program_options.hpp>
#include <fstream>
#include <array>
#include "header.h"
#include <boost/algorithm/string/case_conv.hpp>

using namespace std;
static unsigned g_kmer_len = 5;
static uint8_t g_code[256];
static Embedding g_e;
uint32_t g_mask;

enum class IType {
	Bonito,
	Guppy
};
IType g_itype = IType::Bonito;

void make_code(void) {
	for (int i = 0; i < 256; i++)
		g_code[i] = 4;

	g_code['A'] = g_code['a'] = 0;
	g_code['N'] = g_code['n'] = 0;
	g_code['C'] = g_code['c'] = 1;
	g_code['G'] = g_code['g'] = 2;
	g_code['T'] = g_code['t'] = 3;
}

struct Read {
	string name, seq, qual;
	map<uint32_t, vector<unsigned>> index;
};

struct Motif {
	string name, seq;
	vector<string> eseq;
	unsigned pos, edist;

	bool operator()(const Motif &X, const Motif &Y) const {
		return X.pos < Y.pos;
	}
};

istream &operator>>(istream &in, Read &r) {
	if (g_itype == IType::Guppy) {
		string tmp;
		if (!getline(in, r.name)) return in;
		if (!getline(in, r.seq)) return in;
		if (!getline(in, tmp)) return in;
		if (!getline(in, r.qual)) return in;
	} else {
		if (!getline(in, r.name)) return in;
		string seq;
		r.seq = "";
		do {
			getline(in, seq);
			r.seq += seq;
		} while (in && in.peek() != '>');
	}

	return in;
}

istream &operator>>(istream &in, Motif &m) {
	if (!getline(in, m.name)) return in;
	return getline(in >> std::ws, m.seq); //remove the whitespace in the beginning of line
}

void index_read(Read &r) {
	unsigned len = r.seq.size();
	uint32_t k = 0;
	for (unsigned i = 0; i < g_kmer_len - 1; i++) {
		k = (k << 2) + *(g_code + r.seq[i]);
	}

	for (unsigned i = g_kmer_len - 1; i < len; i++) {
		k = (k << 2) + *(g_code + r.seq[i]);
		unsigned key = k & g_mask;
		r.index[key].push_back(i - (g_kmer_len - 1));
	}
}

void seed(const string &motif, Read &r, vector<unsigned> &candidates,
		  unsigned threshold) {
	unsigned len = motif.size();
	vector<unsigned> all_candidates;
	for (unsigned i = 0; i + g_kmer_len <= len; i += g_kmer_len) {
		uint32_t k = 0;
		for (unsigned j = i; j < i + g_kmer_len; j++)
			k = (k << 2) + *(g_code + motif[j]);
		assert(k == (k & g_mask));
		for (auto c = r.index[k].begin(); c != r.index[k].end(); ++c)
			all_candidates.push_back(*c > i ? *c - i : 0);
	}

	sort(all_candidates.begin(), all_candidates.end());
	unsigned ncandidates = all_candidates.size();
	for (unsigned i = 0; i < ncandidates;) {
		unsigned j = i + 1;
		while (j < ncandidates &&
			all_candidates[i] == all_candidates[j])
			++j;
		if (j - i >= threshold)
			candidates.push_back(all_candidates[i]);
		i = j;
	}
}

struct motif_match {
	unsigned pos, edist;
};

motif_match get_best_match(const string &motif, const vector<string> &evec,
						   Read &r, vector<unsigned> candidates) {
	unsigned min_pos = numeric_limits<unsigned>::min();
	unsigned min_edist = evec[0].size();
	for (unsigned pos : candidates) {
		string_view cview(r.seq.data() + pos, motif.size());
		unsigned edist = g_e.embed_compare(cview, evec, min_edist);
		if (edist < min_edist) {
			min_pos = pos;
			min_edist = edist;
		}
	}

	return motif_match{min_pos, min_edist};
}

struct MapStats {
	bool all_mapped;
	unsigned redist;
};

void map_motifs(vector<Motif> &motifs, Read &r, vector<Motif> &decode_motifs) {

	for (unsigned i = 0; i < motifs.size(); i++) {
		vector<unsigned> candidates;
		seed(motifs[i].seq, r, candidates, 2);
		if (candidates.empty())
			seed(motifs[i].seq, r, candidates, 1);

		if (!candidates.empty()) {
			auto[pos, edist] = get_best_match(motifs[i].seq, motifs[i].eseq,
											  r, candidates);
			motifs[i].pos = pos;
			motifs[i].edist = edist;

			Motif &motif = motifs[i];
			decode_motifs.push_back(motif);
		}

	}
}

MapStats map_motifs(vector<Motif> &motifs, Read &r) {
	unsigned redist = 0;
	bool all_mapped = true;
	for (unsigned i = 0; i < motifs.size(); i++) {
		vector<unsigned> candidates;
		seed(motifs[i].seq, r, candidates, 2);
		if (candidates.empty())
			seed(motifs[i].seq, r, candidates, 1);

		if (!candidates.empty()) {
			auto[pos, edist] = get_best_match(motifs[i].seq, motifs[i].eseq,
											  r, candidates);
			motifs[i].pos = pos;
			motifs[i].edist = edist;
		} else {
			motifs[i].pos = numeric_limits<unsigned>::max();
			motifs[i].edist = motifs[i].eseq.size();
			all_mapped = false;
		}

		redist += motifs[i].edist;
	}

	if (!all_mapped)
		redist = motifs[0].eseq[0].size() * motifs.size();

	return MapStats{all_mapped, redist};
}


void print_decoded_read(vector<Motif> &motifs, Read &r, ostream &out) {

	out << r.name << "\t";
	for (Motif &m : motifs) {
		out << m.name << "," << m.pos << "-";
	}
	out << endl;
}

void print_motifs(vector<Motif> &motifs, ostream &out) {
	sort(motifs.begin(), motifs.end(),
		 [](Motif &left, Motif &right) { return left.pos < right.pos; });

	for (Motif &m : motifs) {
		out << m.name << ",";
		if (m.pos == (numeric_limits<unsigned>::max()))
			out << "*,*\n";
		else
			out << m.pos << "," << m.edist << endl;
	}
}

void print_per_read_motifs(vector<Motif> &motifs, Read &r, ostream &out) {
	sort(motifs.begin(), motifs.end(),
		 [](Motif &left, Motif &right) { return left.pos < right.pos; });

	out << "Read: " << r.name << endl;
	for (Motif &m : motifs) {
		out << m.name << ",";
		if (m.pos == (numeric_limits<unsigned>::max()))
			out << "*,*\n";
		else
			out << m.pos << "," << m.edist << endl;
	}
}

struct ReadStats {
	int ntotal, nmapped;
};

ReadStats decode_reads(vector<Motif> &mmotifs, const string &rfname, ostream &out) {
	ifstream in(rfname);
	if (!in.is_open()) {
		cerr << "ERROR: Invalid input file " << rfname << endl;
		return ReadStats{0, 0};
	}
	cerr << "Processing read file " << rfname << endl;

	vector<Read> reads;
	int ntotal = 0, nmapped = 0;
	do {
		Read r{};
		in >> r;
		ntotal++;

		//read must be atleast as long as a motif
		if (r.seq.size() <= mmotifs[0].seq.size())
			continue;

		index_read(r);

		// find best position for map motifs
		vector<Motif> decode_motifs;
		map_motifs(mmotifs, r, decode_motifs);
		sort(decode_motifs.begin(), decode_motifs.end(), Motif());

		print_decoded_read(decode_motifs, r, out);
	} while (in);

	return ReadStats{ntotal, nmapped};
}

ReadStats process_reads(vector<Motif> &mmotifs, vector<Motif> &qmotifs,
						vector<Motif> &best_mmotifs, Motif &best_qmotif,
						const string &rfname, ostream &out) {
	ifstream in(rfname);
	if (!in.is_open()) {
		cerr << "ERROR: Invalid input file " << rfname << endl;
		return ReadStats{0, 0};
	}
	cerr << "Processing read file " << rfname << endl;

	vector<Read> reads;
	unsigned best_qredist = mmotifs.size() * mmotifs[0].eseq[0].size();
	unsigned best_mredist = best_qredist;
	unsigned best_redist = best_mredist + best_qredist;
	int ntotal = 0, nmapped = 0;
	do {
		Read r{};
		in >> r;
		ntotal++;

		//read must be atleast as long as a motif
		if (r.seq.size() <= mmotifs[0].seq.size())
			continue;

		index_read(r);

		if (qmotifs.empty()) {
			// find best position for map motifs
			auto[all_mapped, mredist] = map_motifs(mmotifs, r);
			if (all_mapped)
				nmapped++;

			if (mredist < best_mredist) {
				best_mredist = mredist;
				best_mmotifs.clear();
				best_mmotifs.insert(best_mmotifs.end(), mmotifs.begin(), mmotifs.end());
			}
		} else {
			// go over qmotif list and find best query motif for this read
			vector<Motif> vq;
			unsigned best_qmotif_idx = numeric_limits<unsigned>::max();
			best_qredist = mmotifs.size() * mmotifs[0].eseq[0].size();
			unsigned best_qpos = 0;
			for (unsigned i = 0; i < qmotifs.size(); ++i) {
				Motif &q = qmotifs[i];
				vq.push_back(q);
				auto[all_mapped, qredist] = map_motifs(vq, r);
				if (all_mapped && qredist < best_qredist) {
					best_qredist = qredist;
					best_qmotif_idx = i;
					best_qpos = vq[0].pos;
					/*
					cerr << "saving qmotif at pos " << best_qpos <<
					" with edist " << best_qredist <<
					" for mredist " << best_mredist <<
					" and total edist " << best_mredist + best_qredist <<
					" when best total was " << best_redist << endl;
					*/
				}
				vq.pop_back();
			}
			if (best_qmotif_idx == numeric_limits<unsigned>::max())
				continue;

			// find best position for map motifs
			auto[all_mapped, best_mredist] = map_motifs(mmotifs, r);
			if (all_mapped)
				nmapped++;

			// if combination of best + best is best, keep both
			if (best_redist > best_mredist + best_qredist) {
				best_redist = best_mredist + best_qredist;

				// save map motifs
				best_mmotifs.clear();
				best_mmotifs.insert(best_mmotifs.end(), mmotifs.begin(), mmotifs.end());

				// save query motifs
				best_qmotif = qmotifs[best_qmotif_idx];
				best_qmotif.pos = best_qpos;
				best_qmotif.edist = best_qredist;

				//cerr << "Finally saving qmotif at pos " << best_qpos <<
				//    " with edist " << best_qredist << endl;
			}
		}
//		print_per_read_motifs(motifs, r, out);
	} while (in);

	return ReadStats{ntotal, nmapped};
}

int main(int argc, char *argv[]) {
	namespace po = boost::program_options;
	bool help{};
	po::options_description description{"motif-search [options]"};
	description.add_options()
		("help,h", po::bool_switch(&help), "Display help")
		("motifs,m", po::value<string>(), "File containing motifs to be mapped")
		("query,q", po::value<string>(), "File containing query motifs")
		("read,r", po::value<vector<string>>()->multitoken(), "File(s) containing reads")
		("kmerlen,l", po::value<unsigned>(), "Length of kmer (5)")
		("output,o", po::value<string>(), "Ouptut file(default stdout)")
		("itype,i", po::value<string>(), "input format (bonito or guppy bc)");
	po::command_line_parser parser{argc, argv};
	parser.options(description);
	auto parsed_result = parser.run();
	po::variables_map vm;
	po::store(parsed_result, vm);
	po::notify(vm);

	if (help) {
		cerr << description << endl;
		return 0;
	}

	if (vm["motifs"].empty() || vm["read"].empty()) {
		cout << "Must specificy motif and atleast one read file as input" << endl;
		return 1;
	}

	if (!vm["kmerlen"].empty()) {
		g_kmer_len = vm["kmerlen"].as<unsigned>();
		cerr << "Using kmer of size " << g_kmer_len << endl;
	}

	assert(g_kmer_len < 16);
	g_mask = (1U << (g_kmer_len * 2)) - 1;

	const string &motif_name = vm["motifs"].as<string>();
	cerr << "Reading map motif file " << motif_name << endl;
	ifstream mfile(motif_name);
	vector<Motif> mmotifs;
	Motif m;
	while (mfile >> m)
		mmotifs.push_back(m);
	cerr << "Found " << mmotifs.size() << " motifs." << endl;

	vector<Motif> qmotifs;
	if (!vm["query"].empty()) {
		const string &query_name = vm["query"].as<string>();
		cerr << "Reading query motif file " << query_name << endl;
		ifstream qfile(query_name);
		while (qfile >> m)
			qmotifs.push_back(m);
		cerr << "Found " << qmotifs.size() << " query motifs." << endl;
	} else {
		cerr << "No query motifs supplied.\n";
	}

	for (Motif &m : mmotifs)
		g_e.embed_string(m.seq, m.eseq);
	for (Motif &q : qmotifs)
		g_e.embed_string(q.seq, q.eseq);

	make_code();

	string ofile_name = "";
	ofstream ofile;
	if (!vm["output"].empty()) {
		ofile_name = vm["output"].as<string>();
		ofile.open(ofile_name);
	}

	if (!vm["itype"].empty()) {
		if (boost::algorithm::to_lower_copy(vm["itype"].as<string>()) ==
			"bonito") {
			cerr << "Setting input format to bonito\n";
			g_itype = IType::Bonito;
		}
	}

	vector<Motif> best_mmotifs;
	Motif best_qmotif;
	auto start = chrono::system_clock::now();
	vector<tuple<string, int, int>> stats;
	for (const string &rf : vm["read"].as<vector<string>>()) {
		auto[ntotal, nmapped] = decode_reads(mmotifs, rf, ofile.is_open() ? ofile : cout);
//		auto[ntotal, nmapped] = process_reads(mmotifs, qmotifs,
//											  best_mmotifs, best_qmotif,
//											  rf, ofile.is_open() ? ofile : cout);
		stats.push_back(tuple{rf, ntotal, nmapped});
	}
	if (!qmotifs.empty())
		best_mmotifs.push_back(best_qmotif);

	if (ofile.is_open())
		ofile.close();

	cerr << "-------------\n";
	cerr << "Summary Stats\n";
	cerr << "-------------\n";
	for (auto &t : stats) {
		cerr << get<2>(t) << " out of " << get<1>(t) << " reads fully mapped in"
														" file " << get<0>(t) << endl;
	}

	auto end = chrono::system_clock::now();
	auto elapsed = chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cerr << "Completed. Wall clock time: " <<
		 elapsed.count() << " millisecs.\n";

	return 0;
}
