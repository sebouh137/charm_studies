

OR=3.1/2
thickness=0.05
z_dome=-2.5
z_tube=-3.75
tube_length=2.5
print(f"""DomeOuter | root | Outer hemisphere shell | 0*cm 0*cm {z_dome}*cm | 0*deg 0*deg 0*deg | aaaaaa | Sphere | {OR-thickness}*cm {OR}*cm 0*deg 360*deg 0*deg 90*deg | G4_Al | no | 1 | 1 | 1 | 1 | 1 | no | no | no
BeampipeTube | root | Hollow pipe | 0*cm 0*cm {z_tube}*cm | 0*deg 0*deg 0*deg | aaaaaa | Tube | {OR-thickness}*cm {OR}*cm {tube_length/2}*cm 0*deg 360*deg | G4_Al | no | 1 | 1 | 1 | 1 | 1 | no | no | no
""")
