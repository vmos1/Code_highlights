L=32
blk=2
m=-0.02
for lev in 0 1 2 3 4;
    do echo $L,3, $blk, 0, $m,$lev;
    ../code/./a1 $L 3 $blk 0 $m $lev 0 0;
    mv results_phi.txt results_phi_L"$L"_m"$m"_lvls"$lev".txt
    mv results_res_lvl-0.txt results_residue_L"$L"_m"$m"_nlvls"$lev".txt
done
