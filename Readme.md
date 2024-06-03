###
Create and run the docker file:
```
docker build --no-cache -t ext .

docker run -it ext
```
Run the experiments:

```
source activate quark_install

cd reduc/

python instance.py 7 7 False True
```

Create Plots:
```
Rscript createPaperPlots.r
```

Load the docker image
```
docker load -i docker_img
```


