# Building Libra's documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/).
Install the pinned documentation dependencies:


```bash
python -m pip install -r requirements.txt
```


Build the static website:
```bash
make html
```

Open `build/html/index.html` in a browser. The build regenerates the complete
`libra_py` API catalog from the source tree. If C++ bindings changed, rebuild
Libra first so the `liblibra_core` page reflects the current Boost.Python
extension.

For a clean warning audit, run:

```bash
make clean html SPHINXOPTS="-W --keep-going"
```
