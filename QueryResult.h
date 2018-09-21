#pragma once
#include <vector>

typedef std::vector<bool> CPURangeQueryResult;

typedef struct RangeQueryResultGPU {
	short jobID;
	short idx;
}RangeQueryResultGPU;

class QueryResult
{
public:
	QueryResult();
	~QueryResult();

	CPURangeQueryResult* start;
};

