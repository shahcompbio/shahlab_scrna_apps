REGISTRY=$1
ORG=$2
VERSION=$3

echo "\n PREPARING \n"

echo "VERSION: ${VERSION}"

echo "\n BUILDING \n"

docker build -t $REGISTRY/$ORG/gex_analysis:$VERSION ../../apps/gex_analysis --no-cache
