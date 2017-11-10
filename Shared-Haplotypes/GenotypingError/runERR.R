#
# analyse ESH errors
#


index.results = function(res.dir, prefix) {
	files = dir(path = res.dir, pattern = sprintf("^%s\\..+\\.RData$", prefix))
	index = data.frame(file = files, 
										 type = sub(sprintf("^%s\\.([^.]+)\\.([0-9\\.]+)\\.([^.]+)\\.RData$", prefix), "\\1", basename(files)),
										 rate = as.numeric(sub(sprintf("^%s\\.([^.]+)\\.([0-9\\.]+)\\.([^.]+)\\.RData$", prefix), "\\2", basename(files))) / 100,
										 indv = sub(sprintf("^%s\\.([^.]+)\\.([0-9\\.]+)\\.([^.]+)\\.RData$", prefix), "\\3", basename(files)),
										 stringsAsFactors = FALSE)
	index
}













