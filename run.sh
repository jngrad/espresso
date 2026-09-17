gmx grompp -f minim.mdp -c polymer.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em


gmx grompp -f spe.mdp -c polymer.gro -p topol.top -o spe.tpr
gmx mdrun -s spe.tpr -o spe.trr -e spe.edr -g spe.log
gmx energy -f spe.edr -s spe.tpr -o pressure_tensor.xvg
20 21 22 23 24 25 26 27 28


22 23 24 25 26 27 28 29 30

gmx traj -f spe.trr -s spe.tpr -of forces.xvg







gmx convert-tpr -s dihedral.tpr -o dihedral2.tpr -generate_velocities -velocity_temp 0
gmx grompp -f spe.mdp -c dihedral.gro -p dihedral.top -o dihedral.tpr
gmx convert-trj -f dihedral.gro -o dihedral.trr -s dihedral.tpr -timestep 0.01
gmx mdrun -rerun dihedral.trr -s dihedral.tpr -o dihedral-rerun.trr -e dihedral-rerun.edr -g dihedral-rerun.log
gmx energy -f dihedral-rerun.edr -s dihedral.tpr -o dihedral-rerun.xvg
gmx dump -f dihedral.trr | less


gmx grompp -f rerun.mdp -c dihedral.gro -p dihedral.top -o rerun.tpr -maxwarn 1
gmx mdrun -rerun dihedral.trr -s rerun.tpr -o forces.trr -e rerun.edr




# that part here works:
mkdir dihedral
python run.py
cd dihedral
# takes 5 min
time for i in {00000..03199}; do
  gmx grompp -f ../spe.mdp -c frame$i.gro -p ../dihedral.top -o frame$i.tpr
  gmx mdrun -s frame$i.tpr -o frame$i.trr -e frame$i.edr -g frame$i.log
  rm confout.gro
  gmx energy -f frame$i.edr -s frame$i.tpr -o frame$i.xvg << EOF
5 10 20 21 22 23 24 25 26 27 28
EOF
done
echo "        Time   Potential    Pressure     Pres-XX     Pres-XY     Pres-XZ     Pres-YX     Pres-YY     Pres-YZ     Pres-ZX     Pres-ZY     Pres-ZZ" > results.dat
tail -n 1 frame*.xvg >> results.dat
sed -i "/==> /d" results.dat
sed -ri "/^$/d" results.dat
