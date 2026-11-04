const $=id=>document.getElementById(id);
const atomic={H:1.008,Li:6.94,C:12.011,N:14.007,O:15.999,Na:22.990,Mg:24.305,Al:26.982,Si:28.085,P:30.974,S:32.06,Cl:35.45,K:39.098,Ca:40.078,Ti:47.867,V:50.942,Cr:51.996,Mn:54.938,Fe:55.845,Co:58.933,Ni:58.693,Cu:63.546,Zn:65.38,Sr:87.62,Y:88.906,Zr:91.224,Nb:92.906,Ba:137.327,La:138.905,Ce:140.116,W:183.84};
function molarMass(s){if(!s)throw Error("Enter the chemical formula.");let i=0;function number(){let m=s.slice(i).match(/^(?:\d+(?:\.\d*)?|\.\d+)/);if(!m)return 1;i+=m[0].length;return +m[0];}function group(){let total=0;while(i<s.length&&s[i]!==")"){if(s[i]==="("){i++;let v=group();if(s[i++]!==")")throw Error("Unmatched parenthesis in formula.");total+=v*number();}else{let m=s.slice(i).match(/^([A-Z][a-z]?)(\d*(?:\.\d+)?)/);if(!m||!(m[1] in atomic))throw Error("Unsupported element or invalid formula near: "+s.slice(i));i+=m[0].length;total+=atomic[m[1]]*(m[2]?+m[2]:1);}}return total;}let m=group();if(i!==s.length||m<=0)throw Error("Invalid chemical formula.");return m;}
function solve(A,b){let n=b.length,M=A.map((r,i)=>[...r,b[i]]);for(let k=0;k<n;k++){let p=k;for(let i=k+1;i<n;i++)if(Math.abs(M[i][k])>Math.abs(M[p][k]))p=i;if(Math.abs(M[p][k])<1e-14)throw Error("Insufficient independent reflections: least-squares fit is rank deficient.");[M[k],M[p]]=[M[p],M[k]];let q=M[k][k];for(let j=k;j<=n;j++)M[k][j]/=q;for(let i=0;i<n;i++)if(i!==k){q=M[i][k];for(let j=k;j<=n;j++)M[i][j]-=q*M[k][j];}}return M.map(r=>r[n]);}
function inv3(m){let d=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);let minor2=m[0][0]*m[1][1]-m[0][1]*m[1][0];if(m[0][0]<=0||minor2<=0||d<=0)throw Error("Fitted reciprocal metric is not positive definite; inspect indexing, crystal system and d values.");return [[(m[1][1]*m[2][2]-m[1][2]*m[2][1])/d,(m[0][2]*m[2][1]-m[0][1]*m[2][2])/d,(m[0][1]*m[1][2]-m[0][2]*m[1][1])/d],[(m[1][2]*m[2][0]-m[1][0]*m[2][2])/d,(m[0][0]*m[2][2]-m[0][2]*m[2][0])/d,(m[0][2]*m[1][0]-m[0][0]*m[1][2])/d],[(m[1][0]*m[2][1]-m[1][1]*m[2][0])/d,(m[0][1]*m[2][0]-m[0][0]*m[2][1])/d,(m[0][0]*m[1][1]-m[0][1]*m[1][0])/d]];}
function fitCell(rows,system){const terms=rows.map(({h,k,l,d})=>{let y=1/(d*d),x;switch(system){case"cubic":x=[h*h+k*k+l*l];break;case"tetragonal":x=[h*h+k*k,l*l];break;case"orthorhombic":x=[h*h,k*k,l*l];break;case"hexagonal":x=[4*(h*h+h*k+k*k)/3,l*l];break;case"rhombohedral":x=[h*h+k*k+l*l,2*(h*k+h*l+k*l)];break;case"monoclinic":x=[h*h,k*k,l*l,2*h*l];break;default:x=[h*h,k*k,l*l,2*h*k,2*h*l,2*k*l];}return{x,y};});let n=terms[0].x.length;if(rows.length<n)throw Error("This crystal system needs at least "+n+" independent reflections; enter more hkl/d rows.");let A=Array.from({length:n},()=>Array(n).fill(0)),b=Array(n).fill(0);for(const t of terms)for(let i=0;i<n;i++){b[i]+=t.x[i]*t.y;for(let j=0;j<n;j++)A[i][j]+=t.x[i]*t.x[j];}let p=solve(A,b),G;if(system==="cubic")G=[[p[0],0,0],[0,p[0],0],[0,0,p[0]]];else if(system==="tetragonal")G=[[p[0],0,0],[0,p[0],0],[0,0,p[1]]];else if(system==="orthorhombic")G=[[p[0],0,0],[0,p[1],0],[0,0,p[2]]];else if(system==="hexagonal")G=[[4*p[0]/3,2*p[0]/3,0],[2*p[0]/3,4*p[0]/3,0],[0,0,p[1]]];else if(system==="rhombohedral")G=[[p[0],p[1],p[1]],[p[1],p[0],p[1]],[p[1],p[1],p[0]]];else if(system==="monoclinic")G=[[p[0],0,p[3]],[0,p[1],0],[p[3],0,p[2]]];else G=[[p[0],p[3],p[4]],[p[3],p[1],p[5]],[p[4],p[5],p[2]]];let D=inv3(G),a=Math.sqrt(D[0][0]),b0=Math.sqrt(D[1][1]),c=Math.sqrt(D[2][2]),deg=x=>Math.acos(Math.max(-1,Math.min(1,x)))*180/Math.PI,alpha=deg(D[1][2]/(b0*c)),beta=deg(D[0][2]/(a*c)),gamma=deg(D[0][1]/(a*b0)),V=Math.sqrt(D[0][0]*(D[1][1]*D[2][2]-D[1][2]*D[2][1])-D[0][1]*(D[1][0]*D[2][2]-D[1][2]*D[2][0])+D[0][2]*(D[1][0]*D[2][1]-D[1][1]*D[2][0]));if(![a,b0,c,V].every(Number.isFinite)||V<=0)throw Error("The fitted cell is physically invalid.");let residuals=terms.map(t=>{let calc=t.x.reduce((s,x,i)=>s+x*p[i],0),dcalc=1/Math.sqrt(calc);return {obs:1/Math.sqrt(t.y),calc:dcalc,res:dcalc-1/Math.sqrt(t.y)};});let rms=Math.sqrt(residuals.reduce((s,r)=>s+r.res*r.res,0)/rows.length);return {a,b:b0,c,alpha,beta,gamma,V,rms,residuals};}
let rows=[{h:1,k:0,l:0,d:4.012},{h:1,k:1,l:0,d:2.83691241},{h:1,k:1,l:1,d:2.31632928},{h:2,k:0,l:0,d:2.006}];let cell=null,rhoTheory=null;
function renderRows(){ $("rows").innerHTML=rows.map((r,i)=>'<div class="hkl-row"><input aria-label="h" data-i="'+i+'" data-k="h" type="number" step="1" value="'+r.h+'"><input aria-label="k" data-i="'+i+'" data-k="k" type="number" step="1" value="'+r.k+'"><input aria-label="l" data-i="'+i+'" data-k="l" type="number" step="1" value="'+r.l+'"><input aria-label="d spacing in angstrom" data-i="'+i+'" data-k="d" type="number" step="any" value="'+r.d+'"><button class="remove" data-remove="'+i+'" aria-label="Remove row">×</button></div>').join("");}

function normalizeHeader(v){return String(v??"").trim().toLowerCase().replace(/[ _-]+/g,"").replace(/[()[\]]/g,"");}
function headerColumn(v){const h=normalizeHeader(v);if(["h","hindex","millerh"].includes(h))return"h";if(["k","kindex","millerk"].includes(h))return"k";if(["l","lindex","millerl"].includes(h))return"l";if(["d","dspacing","dspace","dangstrom","då","d(å)","dvalue","interplanarspacing","spacing"].includes(h))return"d";return null;}
function parseReflectionMatrix(matrix){
  const clean=matrix.filter(row=>row&&row.some(v=>String(v??"").trim()!==""));
  if(!clean.length)throw Error("No reflection data were found in the selected sheet.");
  let headerRow=-1,map={};
  for(let i=0;i<Math.min(clean.length,10);i++){const candidate={};clean[i].forEach((v,j)=>{const key=headerColumn(v);if(key&&candidate[key]===undefined)candidate[key]=j;});if(["h","k","l","d"].every(k=>candidate[k]!==undefined)){headerRow=i;map=candidate;break;}}
  if(headerRow<0){map={h:0,k:1,l:2,d:3};headerRow=-1;}
  const data=clean.slice(headerRow+1);
  const parsed=[];
  for(let ri=0;ri<data.length;ri++){
    const row=data[ri],vals={};
    for(const key of ["h","k","l","d"]){const raw=row[map[key]];const value=typeof raw==="number"?raw:Number(String(raw??"").trim().replace(/,/g,""));if(String(raw??"").trim()===""||!Number.isFinite(value)){if(row.every(v=>String(v??"").trim()===""))continue;throw Error("Invalid or missing "+key+" value in imported data row "+(ri+headerRow+2)+".");}vals[key]=value;}
    if(Object.keys(vals).length)parsed.push(vals);
  }
  if(!parsed.length)throw Error("No valid h, k, l, d-spacing rows were found.");
  if(parsed.some(r=>![r.h,r.k,r.l].every(Number.isInteger)||r.h===0&&r.k===0&&r.l===0||r.d<=0))throw Error("Imported data must contain integer hkl (not all zero) and positive d-spacing values in Å.");
  return parsed;
}
function applyImportedRows(parsed,mode="replace"){
  rows=mode==="append"?rows.concat(parsed):parsed;
  renderRows();
  const status=$("import-status");if(status)status.textContent=parsed.length+" reflection(s) loaded. Review the rows, then fit lattice parameters.";
}
$("rows").addEventListener("paste",e=>{
  const target=e.target;if(!target.matches("input[data-i][data-k]"))return;
  const raw=e.clipboardData?.getData("text");if(!raw||(!raw.includes("\t")&&!raw.includes("\n")))return;
  e.preventDefault();
  try{const matrix=raw.replace(/\r/g,"").split("\n").filter(x=>x.trim()!=="").map(x=>x.split("\t"));const parsed=parseReflectionMatrix(matrix);applyImportedRows(parsed,"replace");}
  catch(err){const status=$("import-status");if(status)status.textContent="Paste error: "+err.message;}
});
$("excel-file").addEventListener("change",async e=>{
  const file=e.target.files?.[0];if(!file)return;
  const status=$("import-status");status.textContent="Reading "+file.name+"…";
  try{
    const ext=file.name.split(".").pop().toLowerCase();let matrix;
    if(ext==="csv"||ext==="tsv"){
      const text=await file.text();matrix=text.replace(/^\uFEFF/,"").split(/\r?\n/).filter(x=>x.trim()!=="").map(line=>{const delimiter=ext==="tsv"||line.includes("\t")?"\t":",";return line.split(delimiter).map(v=>v.trim().replace(/^"(.*)"$/,"$1"));});
    }else{
      if(!window.XLSX)throw Error("Excel reader did not load. Check your internet connection and reload the page.");
      const buffer=await file.arrayBuffer(),workbook=XLSX.read(buffer,{type:"array",cellDates:false});
      const sheetName=workbook.SheetNames[0];if(!sheetName)throw Error("The workbook contains no worksheets.");
      matrix=XLSX.utils.sheet_to_json(workbook.Sheets[sheetName],{header:1,raw:true,defval:""});
      status.textContent="Reading worksheet: "+sheetName+"…";
    }
    const parsed=parseReflectionMatrix(matrix);applyImportedRows(parsed,"replace");
  }catch(err){status.textContent="Import error: "+err.message;}
  e.target.value="";
});
$("rows").addEventListener("input",e=>{let i=e.target.dataset.i,k=e.target.dataset.k;if(i!==undefined)rows[+i][k]=Number(e.target.value);});
$("rows").addEventListener("click",e=>{if(e.target.dataset.remove!==undefined){rows.splice(+e.target.dataset.remove,1);renderRows();}});
$("add-row").onclick=()=>{rows.push({h:0,k:0,l:1,d:1});renderRows();};
$("fit").onclick=()=>{cell=null;rhoTheory=null;$("theory-output").textContent="Fit lattice parameters first.";$("relative-output").textContent="Requires theoretical density and sample mass/volume.";$("bulk-output").textContent="Enter sample mass and volume.";try{if(rows.some(r=>![r.h,r.k,r.l].every(Number.isInteger)||r.h===0&&r.k===0&&r.l===0||!Number.isFinite(r.d)||r.d<=0))throw Error("Every row needs integer hkl (not all zero) and positive d-spacing.");cell=fitCell(rows,$("system").value);$("cell-output").innerHTML='<div class="metrics"><div><small>a</small><b>'+cell.a.toFixed(6)+' Å</b></div><div><small>b</small><b>'+cell.b.toFixed(6)+' Å</b></div><div><small>c</small><b>'+cell.c.toFixed(6)+' Å</b></div><div><small>α</small><b>'+cell.alpha.toFixed(5)+'°</b></div><div><small>β</small><b>'+cell.beta.toFixed(5)+'°</b></div><div><small>γ</small><b>'+cell.gamma.toFixed(5)+'°</b></div><div><small>Unit-cell volume</small><b>'+cell.V.toFixed(5)+' Å³</b></div><div><small>d-spacing RMS</small><b>'+cell.rms.toExponential(3)+' Å</b></div></div><p class="hint">Least-squares fit of 1/d² against the reciprocal metric for the selected crystal system.</p><table><thead><tr><th>Reflection</th><th>d observed (Å)</th><th>d calculated (Å)</th><th>Residual (Å)</th></tr></thead><tbody>'+rows.map((r,i)=>'<tr><td>'+r.h+r.k+r.l+'</td><td>'+r.d.toFixed(5)+'</td><td>'+cell.residuals[i].calc.toFixed(5)+'</td><td>'+cell.residuals[i].res.toExponential(3)+'</td></tr>').join("")+'</tbody></table>';updateDensity();}catch(e){cell=null;rhoTheory=null;$("cell-output").textContent=e.message;$("theory-output").textContent="Fit valid lattice parameters first.";$("relative-output").textContent="Requires valid theoretical density and sample mass/volume.";$("bulk-output").textContent="Enter valid sample mass and volume.";}};
function updateDensity(){if(!cell){rhoTheory=null;$("theory-output").textContent="Fit lattice parameters first.";$("relative-output").textContent="Requires theoretical density and sample mass/volume.";return;}try{let M=molarMass($("formula").value),Z=Number($("z").value);if(!Number.isFinite(Z)||Z<=0)throw Error("Z must be a finite positive number.");rhoTheory=Z*M/(6.02214076e23*cell.V*1e-24);if(!Number.isFinite(rhoTheory)||rhoTheory<=0)throw Error("Theoretical density is not finite; check formula, Z and cell volume.");$("theory-output").innerHTML='<div class="big-result">'+rhoTheory.toFixed(6)+' <small>g cm⁻³</small></div><p>Molar mass: '+M.toFixed(5)+' g mol⁻¹ · Z = '+Z+' · V = '+cell.V.toFixed(5)+' Å³</p>';updateBulkAndRelative();}catch(e){rhoTheory=null;$("theory-output").textContent=e.message;$("relative-output").textContent="Requires a valid theoretical density and sample mass/volume.";}}
$("formula").addEventListener("input",updateDensity);$("z").addEventListener("input",updateDensity);
$("calc-theory").onclick=updateDensity;
function updateBulkAndRelative(){let m=Number($("mass").value),V=Number($("bulk-volume").value);if(Number.isFinite(m)&&Number.isFinite(V)&&m>0&&V>0){let bulk=m/V;$("bulk-output").innerHTML='<div class="big-result">'+bulk.toFixed(6)+' <small>g cm⁻³</small></div>';if(Number.isFinite(rhoTheory)&&rhoTheory>0){let rel=100*bulk/rhoTheory;$("relative-output").innerHTML='<div class="metrics"><div><small>Bulk density</small><b>'+bulk.toFixed(5)+' g cm⁻³</b></div><div><small>Relative density</small><b>'+rel.toFixed(3)+'%</b></div><div><small>Estimated porosity</small><b>'+(100-rel).toFixed(3)+'%</b></div></div>'+(rel>100?'<p class="warn">Relative density exceeds 100%. Check the phase composition, Z, XRD fit, mass and volume.</p>':'');}else{$("relative-output").textContent="Requires a valid theoretical density.";}}else{$("bulk-output").textContent="Enter positive mass (g) and volume (cm³).";$("relative-output").textContent="Requires valid sample mass/volume and theoretical density.";}}
$("calc-bulk").onclick=updateBulkAndRelative;
$("mass").addEventListener("input",updateBulkAndRelative);$("bulk-volume").addEventListener("input",updateBulkAndRelative);
renderRows();