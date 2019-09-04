#!/usr/bin/awk

BEGIN{}
{
    if ( match($0, /CB:Z:[ACGT]+/) )
    {
        REF = $4 - 1;
        QUE = 0;
        CIG = $6;
        CEL = substr($0, RSTART+5, RLENGTH-5);
        while( match(CIG, /^[[:digit:]]+/) )
        {
            N = substr(CIG, RSTART, RLENGTH);
            CIG = substr(CIG, RSTART+RLENGTH);
            if( match(CIG, /^[MIDNSHP=X]/) )
            {
                C = substr(CIG, RSTART, RLENGTH);
                CIG = substr(CIG, RSTART+RLENGTH);
                if (C == "M" || C == "=" || C == "X")
                {                    
                    REF += N;
                    QUE += N;
                    if (TAG <= REF)
                    {
                        X[CEL, substr($10, QUE - REF + TAG, 1)]++;
                        next;
                    };
                } else if (C == "D" || C == "N")
                {
                    REF += N;
                    if (TAG <= REF)
                    {
                        X[CEL, "N"]++;
                        next;
                    }
                } else if (C == "I" || C == "S")
                {
                    QUE += N;
                };
            };
        };
    };
}
END{ for (p in X) { split(p, x, SUBSEP); print x[1], x[2], X[x[1], x[2]] } }

