#pragma once
class Embedding {
	public:
	Embedding();
	void embed_string(const std::string_view &str,
					  std::vector<std::string> &evec);
	unsigned embed_compare(const std::string_view &to_embed,
						   const std::vector<std::string> &evec, unsigned threshold);

	private:
	std::string embed_string(const std::string_view &str, int str_id);
	unsigned embed_compare(const std::string_view &to_embed,
						   const std::string &ref, int str_id, unsigned threshold);

	constexpr static int N_RNDSTR = 8;
	constexpr static int N_ECHAR = 5;
	constexpr static int MAX_ELEN = 100;
	constexpr static int N_RBITS_PER_STR = MAX_ELEN * N_ECHAR;
	constexpr static int N_RNDBITS = N_RBITS_PER_STR * N_RNDSTR;
	constexpr static int EFACTOR = 3;
	constexpr static int EPAD = 4;
	std::bitset<N_RNDBITS> rnd_str;
};

#define BITPOS(STR_ID, OFFSET, CHAR_ID) (STR_ID * N_RBITS_PER_STR +\
        OFFSET * N_ECHAR + CHAR_ID)

