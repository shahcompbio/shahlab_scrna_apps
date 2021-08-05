
# Prereq: build.sh run

REGISTRY=$1
ORG=$2
VERSION=$3

echo "\n PREPARING TO TEST \n"

docker exec $REGISTRY/$ORG/gex_analysis:$VERSION pytest --pyargs gex_analysis
