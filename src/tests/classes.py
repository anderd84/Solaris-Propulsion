from Cooling import domain
from Cooling import material

material.DomainMaterial

print(material.DomainMaterial.FREE == 0)

coolmesh: domain.DomainMC = domain.DomainMC.LoadFile("coolmesh.msh")

test = domain.DomainMMAP(coolmesh)

print(test.temperature[0,0])
test.x[0,0] = 1
print(test.x[0,0])
print(test.previousFlow[0,0])
test.previousFlow[0,0] = [1,2]
print(test.previousFlow[0,0])

test.material[0,0] = material.DomainMaterial.COWL
print(test.material[0,0])