#!/usr/bin/gawk -f

# INFO:Matrix_A:info:           1      104286       35937     1307637      134532
# INFO:Matrix_A:Rows:            1           1           1
# INFO:Matrix_A:Cols_Value:            1           1  0.28609781383975469
# INFO:Vector_b:            1           1   5.5114984137158939E-004

# INFO:Matrix_A:info:           1      176862      176862     4787546      176862
# INFO:Matrix_A:Rows:            0           0
# INFO:Matrix_A:Cols_Value:            0           0   1.00000000000000     
# INFO:Vector_b:Values:            0  0.000000000000000E+000

 
/^ INFO:Matrix_A:info:/ {
    printf("# %d %d %d\n", $2, $3, $4);
}

#/^ INFO:Number/ {
#    printf("# %d %d %d\n", idx0, idx1, $8);
#}

/^ INFO:Matrix_A:Rows:/ {
    print "A_rows", $2, $3
}

/^ INFO:Matrix_A:Cols_Value:/ {
    print "A_cols", $2, $3, $4
}

/^ INFO:Vector_b:/ {
    print "b", $2, $3
}

/^ INFO:Vector_u:/ {
    print "u", $2, $3
}

/^ INFO:Vector_v:/ {
    print "v", $2, $3
}

/^ INFO:Vector_w:/ {
    print "w", $2, $3
}

/^ INFO:Vector_x:/ {
    print "x", $2, $3
}

/^ INFO:Vector_u1:/ {
    print "u1", $2, $3
}

/^ INFO:Vector_u2:/ {
    print "u2", $2, $3
}

/^ INFO:Vector_u3:/ {
    print "u3", $2, $3
}
