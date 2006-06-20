\set semicolon_hack=1
\go

--
-- featureloc fmin-fmax checks
--
select 'Checking featureloc records where fmin>fmax';

select count(*)
from featureloc
where fmin>fmax;


--
-- dangling feature_dbxref records
--
select 'Checking for dangling feature_dbxref records';

select count(fd.feature_id) 
from feature_dbxref fd 
where not exists (
select 1 
from feature f 
where f.feature_id = fd.feature_id )

--
-- dangling organism records
--
select 'Checking for dangling organism records';

select o.genus, o.species, count(o.organism_id) 
from organism o 
where not exists (
select 1 
from feature f 
where f.organism_id = o.organism_id)
group by o.genus, o.species
having count(o.organism_id) > 0
