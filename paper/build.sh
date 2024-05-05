pandoc paper.md --out paper.tex --biblatex
curl -H "Authorization: Bearer $SCIWHEEL_API_TOKEN" 'https://sciwheel.com/extapi/work/references/export?projectId=879256' > paper.bib
cd template_frontiers
latexmk frontiers.tex
