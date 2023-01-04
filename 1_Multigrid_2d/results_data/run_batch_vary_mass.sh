L=32
blk=2
#for lev in 0 2 3 4;
for lev in 1;
    do for m in -0.022 -0.015 -0.01 -0.005;
        do echo $L,3, $blk, 0, $m,$lev;
        ../code/./a1 $L 3 $blk 0 $m $lev 0 0;
        mv results_phi.txt results_phi_L"$L"_m"$m"_lvls"$lev".txt
         mv results_res_lvl-0.txt results_residue_L"$L"_m"$m"_nlvls"$lev".txt
    done
done
