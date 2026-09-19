"""Cover magnetic rectangle using Weyl/Davis-Kahan bounds on projectors.

Analytical enclosure with floating-point eigensolvers and explicit numerical
padding, not an interval-arithmetic/formal proof checker. If the box budget
is exhausted, the retained upper bound remains conservative in exact arithmetic.
"""
import heapq,itertools,json
import numpy as np
import nv_model as nv
from search_interference_temperature import OUT

V=nv.Hes((0,0,0),0)-np.kron(nv.I2,nv.Hgs((0,0,0)))
w=np.linalg.eigvalsh(V);V=V-(w[0]+w[-1])/2*np.eye(6)
VMAX=float((w[-1]-w[0])/2)+1e-10
PERMS=list(itertools.permutations(range(3)))


def enclosure(box):
    ax,bx,az,bz=box
    x=(ax+bx)/2;z=(az+bz)/2
    w,U=np.linalg.eigh(nv.Hgs((x,0,z)))
    eps=nv.GE*np.hypot((bx-ax)/2,(bz-az)/2)+1e-12
    gaps=np.array([min(abs(w[i]-w[j]) for j in range(3) if j!=i) for i in range(3)])
    errs=np.minimum(1.,eps/np.maximum(gaps-eps,1e-300))
    scores=[sum(abs(U[bare,p[bare]])**2 for bare in range(3)) for p in PERMS]
    good=max(s for p,s in zip(PERMS,scores) if p[1]==0)
    other=max(s for p,s in zip(PERMS,scores) if p[1]!=0)
    assignment_certified=bool(good-other>2*sum(errs))
    P=np.kron(nv.I2,np.outer(U[:,0],U[:,0].conj()))
    local=float(np.linalg.norm(P@V@(np.eye(6)-P),2))
    bound=min(VMAX,local+2*VMAX*errs[0]+1e-10) if assignment_certified else VMAX
    return bound,local,assignment_certified


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    root=(0.,.5,.005,.05);u,l,_=enclosure(root)
    heap=[(-u,0,root)];best=l;counter=1;processed=0
    target=4.13
    while heap and -heap[0][0]>target and processed<200000:
        neg,_,box=heapq.heappop(heap);processed+=1
        ax,bx,az,bz=box
        if bx-ax>=bz-az:
            m=(ax+bx)/2; children=[(ax,m,az,bz),(m,bx,az,bz)]
        else:
            m=(az+bz)/2;children=[(ax,bx,az,m),(ax,bx,m,bz)]
        for child in children:
            ub,lb,ok=enclosure(child);best=max(best,lb if ok else 0.)
            heapq.heappush(heap,(-ub,counter,child));counter+=1
        if processed%10000==0:
            print(json.dumps(dict(boxes=processed,upper=-heap[0][0],lower=best)),flush=True)
    upper=-heap[0][0]
    result=dict(domain=dict(Bx=[0,.5],Bz=[.005,.05]),processed_boxes=processed,
        leaf_boxes=len(heap),mixing_norm_upper_GHz=upper,mixing_norm_found_GHz=best,
        angular_mixing_bound=2*np.pi*upper,target_reached=bool(upper<=target),
        enclosure_method='Weyl + Davis-Kahan + explicit Hungarian assignment margin',
        caveat='floating point with padding; not outward-rounded interval arithmetic')
    (OUT/'magnetic_enclosure.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    # Keep the complete covering so every box can be independently checked.
    np.savez_compressed(OUT/'magnetic_cover.npz',boxes=np.array([h[2] for h in heap]),
                        upper=np.array([-h[0] for h in heap]))
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
