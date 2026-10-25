import re

# Read the file
with open('skewImpCorrectedSnGrad.C', 'r') as f:
    content = f.read()

# Fix internal field divisions
content = re.sub(r'/sqr\(magSf\)', '/max(sqr(magSf), SMALL)', content)

# Fix boundary field divisions  
content = re.sub(r'/sqr\(mesh\.boundary\(\)\[patchI\]\.magSf\(\)\)', '/max(sqr(mesh.boundary()[patchI].magSf()), SMALL)', content)

# Write back
with open('skewImpCorrectedSnGrad.C', 'w') as f:
    f.write(content)

print("Fixed division by zero issues")
