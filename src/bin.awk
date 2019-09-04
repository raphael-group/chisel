#!/usr/bin/awk


BEGIN{}
{
    if ( match($0, /CB:Z:[ACGT]+/) )
    {
        X[substr($0, RSTART+5, RLENGTH-5)]++
    }
}
END{ for(i in X) print i, x[i] }
