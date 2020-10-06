


COMPONENTS=$(wildcard src/components/*.js) $(wildcard src/*.html) $(wildcard src/*.js*)
dist=dist/tenx-websummary-script.min.js dist/tenx-websummary-styles.min.css

.PHONY: all
all:$(dist)

$(dist): $(COMPONENTS)
	cd src; yarn install
	cd src; yarn run build

clean:
	rm -rf dist src/node_modules
