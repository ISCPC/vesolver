#!/usr/bin/gawk -f

BEGIN {
    sum = 0;
}

/^TIME:/ {
    print $2;
    sum += $2;
}

END {
    printf("TOTAL: %f [sec]\n", sum);
}
