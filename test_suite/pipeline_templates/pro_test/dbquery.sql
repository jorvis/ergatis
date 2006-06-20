\set semicolon_hack=1
\go
select "Q1: Number of features by type";

select substring(c.name,1,30),count(f.uniquename) 
from feature f, cvterm c 
where f.type_id = c.cvterm_id 
group by c.name order by 1;

select "Q2: Number of features associated with a compute";

select substring(c.name,1,30),substring(a.name,1,30),a.analysis_id,count(f.uniquename) 
from feature f, analysisfeature af, analysis a, cvterm c 
where f.feature_id = af.feature_id 
and af.analysis_id = a.analysis_id 
and f.type_id = c.cvterm_id 
group by c.name,a.name,a.analysis_id order by 1;

select "Q3: Number and feature properties by type";

select substring(fc.name,1,30),substring(c.name,1,30),count(f.uniquename) 
from feature f, featureprop fp, cvterm c, cvterm fc 
where f.feature_id = fp.feature_id 
and fp.type_id = c.cvterm_id 
and f.type_id = fc.cvterm_id 
group by fc.name,c.name order by 1;

select "Q4: Number of features by organism";

select substring(o.common_name,1,30),substring(o.abbreviation,1,30),count(f.uniquename) 
from feature f, organism o 
where f.organism_id = o.organism_id 
group by o.common_name,o.abbreviation order by 1;

select "Q5: Number of dbxref records by external db";

select substring(d.name,1,30),count(x.accession)
from dbxref x, db d
where x.db_id = d.db_id
group by d.name order by 1;

select "Q5: Number of primary database references for all features";

select d.name,count(x.accession)
from dbxref x, db d, feature f
where x.db_id = d.db_id
and f.dbxref_id = x.dbxref_id
group by d.name order by 1;

select "Q6: Number of secondary database references for all features";

select d.name,count(x.accession)
from dbxref x, db d, feature_dbxref f
where x.db_id = d.db_id
and f.dbxref_id = x.dbxref_id
group by d.name order by 1;


select "Q7: F Number of localizations per feature type";

select substring(c.name,1,30),count(f.uniquename)
from feature f, feature fs, featureloc fl, cvterm c
where f.feature_id = fl.feature_id
and fs.feature_id = fl.srcfeature_id
and fs.type_id = c.cvterm_id
group by c.name order by 1;

select "Q8: Number of feature_relationships by type";

select substring(c.name,1,30),count (fg.object_id)
from feature_relationship fg, cvterm c
where c.cvterm_id = fg.type_id
group by c.name order by 1;

select "Q9: Number of features by type by sequence";

select substring(c.name,1,30),fs.uniquename,count(f.uniquename)
from feature f, feature fs, featureloc fl, cvterm c, cvterm a
where f.feature_id = fl.feature_id
and fs.feature_id = fl.srcfeature_id
and f.type_id = c.cvterm_id
and fs.type_id = a.cvterm_id
and a.name = 'assembly'
group by c.name,fs.uniquename order by 1;

select "Q10: Number of secondary feature_cvterms by type";

select substring(c.name,1,30),count(fc.feature_id)
from feature f, feature_cvterm fc, cvterm c, cvterm a
where f.type_id = c.cvterm_id
and fc.feature_id = f.feature_id
and fc.cvterm_id = a.cvterm_id
group by c.name order by 1;

select "Q11: analysisfeature relationships";

select c.name,a.program,count(f.uniquename) 
from feature f, analysisfeature af, cvterm c, analysis a 
where af.type_id = c.cvterm_id 
and f.feature_id = af.feature_id 
and af.analysis_id = a.analysis_id
group by c.name,a.program order by 1;
