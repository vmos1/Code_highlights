L=32
blk=2
m=-0.02
for lev in 2 3 4;
    do for cpy in 1 2 3 4;
        do echo $L,3, $blk, 0, $m,$lev;
        ../code/./a1 $L 3 $blk 0 $m $lev 1 $cpy;
        mv results_phi.txt results_phi_L"$L"_m"$m"_lvls"$lev"_blk"$blk"_cpy"$cpy".txt
        mv results_residue.txt results_residue_L"$L"_m"$m"_lvls"$lev"_blk"$blk"_cpy"$cpy".txt
        done;
        
    ./a1 $L 3 $blk 0 $m $lev 0 0;
    mv results_phi.txt results_phi_L"$L"_m"$m"_lvls"$lev"_blk"$blk"_cpy0.txt
    mv results_residue.txt results_residue_L"$L"_m"$m"_lvls"$lev"_blk"$blk"_cpy0.txt
done
