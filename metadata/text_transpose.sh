# !/bin/bash
# ----------------------------------------------------------------------------
# Copyright (c) 2016--, Gregory Poore.
#
# Distributed under the terms of the Modified BSD License.
#
# Dependencies:
# - GNU datamash (http://www.gnu.org/software/datamash/)
# 
# Resources used in making this script:
# - http://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
# - http://onetipperday.sterding.com/2014/09/transpose-tab-delimited-file.html
# ----------------------------------------------------------------------------

for f in All_CDEs*.txt
do
    echo "Transposing $f file... Saved as Transposed_$f"
    cat $f | datamash transpose > Transposed_$f
done

echo "Transposing complete!"

### Extract headers to another file to find intersection of CDE values.
### After closer examination, it was discovered that not all of the
### metadata tables had the same number of columns. Hence, an intersection
### of the column headers is performed, first by extracting the header names,
### and then by performing set operations in Python.

echo "Extracting column headers and saving in column_names_CDEs.txt ..."

grep '^bcr_patient_barcode' -h Transposed_All_CDEs*.txt > column_names_CDEs.txt

echo "Complete!"