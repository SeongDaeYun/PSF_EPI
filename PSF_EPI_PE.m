
% ========================================================================
%
% This script is written by Seong Dae Yun. 
% 27.06.2022 Seong Dae Yun (s.yun@fz-juelich.de)
%
% * Function: [FWHM, h0] =
%             PSF_EPI_PE( sz_y, T2, TE, espc, pi_f, pf_f, is_T2d, is_PF1 )
%
%             - simulates PSF of EPI in the PE direction
%             - plots |PSF| and calculates its FWHM
% 
%
% * Example :
%   1) FWHM = PSF_EPI_PE( 192, 33.2, 50, 1.2, 3, 6/8 ) ;
%
%
% * Input : ( sz_y, T2, TE, espc, pi_f, pf_f, is_T2d, is_PF1 )
%           --------------------------------------------------------------
%           Input (Default)    |               Description
%           --------------------------------------------------------------
%           sz_y        (N/A)  | Phase encoding size
%           T2          (N/A)  | T2 or T2* time
%           TE          (N/A)  | Echo time (ms)
%           espc        (N/A)  | Echo spacing time (ms)
%           pi_f        (N/A)  | Parallel imaging acceleration factor
%                                e.g.) 1, 2, 3, ...
%           pf_f        (N/A)  | Partial Fourier imaging factor ( >= 0.5 )
%                                e.g.) 5/8, 6/8, ... 
%           is_T2d      (true) | Option to enable/disable T2 decay 
%           is_PF1      (true) | Option to simulate PF trajectory 
%                                true: one-sided PF, false: two-sided PF 
%           --------------------------------------------------------------
%
%
% * Output: [ FWHM, h0 ]
%           --------------------------------------------------------------
%           Output             |               Description
%           --------------------------------------------------------------
%           FWHM               | calculated FWHM (in pixel)
%           h0                 | figure handle 
%           --------------------------------------------------------------
% 
%
% * By using this program, you agree with the following conditions: 
%   1) This script is for non-commercial applications.
%   2) This script is for research purposes, 
%      and should not be used in any diagnostic setting.
%
% ========================================================================




function [FWHM, h0] = PSF_EPI_PE( sz_y, T2, TE, espc, pi_f, pf_f, is_T2d, is_PF1 )


% ========================================================================
% Check input arguments
% ========================================================================
% Check number
if ( nargin < 6 )
    error( 'Not enough input arguments' ) ;
elseif ( nargin == 6 )
    is_T2d = true ;
    is_PF1 = true ;
elseif ( nargin == 7 )
    is_PF1 = true ;
elseif ( nargin > 8 )
    error( 'Too many input arguments' ) ;
end ;

% Check input value
if ( pf_f < 0.5 )
    error( 'Too small partial Fourier factor' ) ;
end ;
if ( ~( is_T2d == true || is_T2d == false ) )
    error( 'is_T2d can be either true or false.' ) ;
end ;
if ( ~( is_PF1 == true || is_PF1 == false ) )
    error( 'is_PF1 can be either true or false.' ) ;
end ;
if ( ceil( sz_y*pf_f/pi_f ) ~= floor( sz_y*pf_f/pi_f ) )
    sz_red = sz_y*pf_f/pi_f ;
    warning( ['The souce size is not completely divided by the acceleration factors: ', num2str( sz_red )] ) ;
end


% ========================================================================
% Define point-source
% ========================================================================
% Image domain
py_ctr = ceil( ( sz_y+1 )/2 ) ;
ps_src = zeros( sz_y, 1 ) ;
ps_src( py_ctr ) = 1 ;

% Frequency doman
ks_src = ifftshift( ifft( ifftshift( ps_src, 1 ), sz_y, 1 ), 1 ) ;
if ( is_PF1 == true )
    ks_idx = ( sz_y - sz_y*pf_f + 1 ):1:sz_y ;
    ky_ctr = py_ctr - ( sz_y - sz_y*pf_f ) ; 
else
    ks_pts = py_ctr - floor( ( sz_y*pf_f ) / 2 ) ;
    ks_idx = ks_pts:1:( ks_pts + sz_y*pf_f - 1 ) ;
    ky_ctr = py_ctr - ks_pts + 1 ;
end ;
len_ks = length( ks_idx ) ;


% ========================================================================
% Apply T2 decay
% ========================================================================
% TE range
espc_e = espc/pi_f ;
TE_i   = 0:espc_e:( len_ks-1 )*espc_e ; 
TE_r   = TE_i + TE - TE_i( ky_ctr ) ;

% Apply
ks_res = zeros( sz_y, 1 ) ;
t2_wgt = ones( len_ks, 1 ) ;
if ( is_T2d )
    t2_wgt = exp( -TE_r/T2 ).' ;
end
ks_res( ks_idx ) = ks_src( ks_idx ).*t2_wgt ;


% ========================================================================
% Result: PSF
% ========================================================================
ps_res  = fftshift( fft( fftshift( ks_res, 1 ), [], 1 ), 1 ) ;
ps_resn = ps_res/max( abs( ps_res ) ) ;


% ========================================================================
% Display
% ========================================================================
[FWHM, x_pos1, x_pos2] = get_fwhm_psf( ps_res ) ;
h0 = figure ;
set( gcf, 'Position', [50 50 1400 1000] ) ;
set( gcf, 'Color', [1 1 1] ) ;

% a) Original point-source
h1 = subplot( 2, 2, 1 ) ; hold on ; set( h1, 'Tag', 'h1' ) ;
plot( 1:1:sz_y, ps_src ) ;
title( 'Point-source' )
xlabel( 'Image pixel index' ) ;
ylabel( 'Magnitude' ) ;
axis( [0 sz_y+1, 0.0 1.1] ) ;

% b) PSF (reconstructed point-source)
h2 = subplot( 2, 2, 2 ) ; hold on ; set( h2, 'Tag', 'h2' ) ;
plot( 1:1:sz_y, abs( ps_resn ) ) ;
plot( [1, sz_y], [0.5, 0.5], 'r--' ) ;
plot( [x_pos1, x_pos2], [0.5, 0.5], 'r.-', 'MarkerSize', 10, 'MarkerEdgeColor','black' ) ;
title( 'Reconstructed point-source (PSF)' )
xlabel( 'Image pixel index' ) ;
ylabel( 'Magnitude' ) ;
axis( [0 sz_y+1, 0 max( abs( ps_resn ) )*1.1] ) ;
legend( 'PSF', 'FWHM line' ) ;

arg1{1} = h2 ;
p1 = [floor( x_pos1-sz_y/8 ) 0.5] ;
p2 = [x_pos1 0.5] ;
[xn, yn] = get_norm_pos( h2, p1, p2 ) ;
hn = annotation( 'Arrow', xn, yn ) ;
set( hn, 'Tag', 'arrow_h2' ) ;
arg1{2} = p1 ; arg1{3} = p2 ;

p1 = [floor( x_pos2+sz_y/8 ) 0.5] ;
p2 = [x_pos2 0.5] ;
[xn, yn] = get_norm_pos( h2, p1, p2 ) ;
hn = annotation( 'Arrow', xn, yn ) ;
set( hn, 'Tag', 'arrow_h2' ) ;
arg1{4} = p1 ; arg1{5} = p2 ;

txt_x = h2.Position(1) + h2.Position(3)*0.05 ;
txt_y = h2.Position(2) + h2.Position(4)*0.50 ;
txt_w = h2.Position(3)*0.25 ;
txt_h = h2.Position(4)*0.08 ;
annotation( 'textbox', [txt_x txt_y txt_w txt_h], 'String', ['FWHM: ', num2str( FWHM )] ) ;

% c) k-space of the point-source
h3 = subplot( 2, 2, 3 ) ; hold on ; set( h3, 'Tag', 'h3' ) ;
plot( 1:1:sz_y, abs( ks_src ) ) ;
title( 'kspace of the point-source' )
xlabel( 'Kspace pixel index' ) ;
ylabel( 'Magnitude' ) ;
axis( [0 sz_y+1, 0 max( abs( ks_src ) )*1.1] ) ;

% d) k-space of the PSF
h4 = subplot( 2, 2, 4 ) ; hold on ; set( h4, 'Tag', 'h4' ) ;
ha = plot( ks_idx, abs( ks_res( ks_idx ) ) ) ;
title( 'kspace of the PSF' )
xlabel( 'Kspace pixel index' ) ;
ylabel( 'Magnitude' ) ;
axis( [0 sz_y+1, 0 max( abs( ks_res ) )*1.1] ) ;

pi = [1 0] ;
p0 = [ks_idx(1)-1 0] ;
p1 = [ks_idx(1) 0] ;
p2 = [ks_idx(end) 0] ;
p3 = [ks_idx(end)+1 0] ;
pe = [sz_y 0] ;
hb = plot( ks_idx(1)     , ks_res( ks_idx(1)   ), 'black>' ) ;
hc = plot( py_ctr        , ks_res( py_ctr      ), 'black*' ) ;
hd = plot( ks_idx(end)   , ks_res( ks_idx(end) ), 'black<' ) ;
he = plot( [p1(1) p2(1)] , [p1(2) p2(2)]        , 'blue'   ) ;

arg1{6} = h4 ;
[xn, yn] = get_norm_pos( h4, p1, p2 ) ;
hn = annotation( 'doublearrow', xn, yn, 'Color', [0 0 1] ) ;
set( hn, 'Tag', 'arrow_h4' ) ;
arg1{7} = p1 ; arg1{8} = p2 ;

if ( pf_f < 1 )
    hf = plot( [pi(1) p0(1)] , [pi(2) p0(2)]    , 'red'    ) ;
    [xn, yn] = get_norm_pos( h4, pi, p0 ) ;
    hn = annotation( 'doublearrow', xn, yn, 'Color', [1 0 0] ) ;
    set( hn, 'Tag', 'arrow_h4' ) ;
    arg1{9} = pi ; arg1{10} = p0 ;
    if ( is_PF1 == false )
        [xn, yn] = get_norm_pos( h4, p3, pe ) ;
        hn = annotation( 'doublearrow', xn, yn, 'Color', [1 0 0] ) ;
        set( hn, 'Tag', 'arrow_h4' ) ;
        arg1{11} = p3 ; arg1{12} = pe ;
    else
        arg1{11} = [0 0] ; arg1{12} = [0 0] ; 
    end 
    
    legend( [ha, hb, hc, hd, hf, he], ...
            'Kspace-signal decay', ...
            ['Echo (', num2str(TE_r(1))     , ') ms @Start'], ...
            ['Echo (', num2str(TE_r(ky_ctr)), ') ms @TE'   ], ...
            ['Echo (', num2str(TE_r(end))   , ') ms @End'  ], ...
            'Skipped by PF', 'Acquisition window' ) ;
else
    legend( [ha, hb, hc, hd, he], ...
            'Kspace-signal decay', ...
            ['Echo (', num2str(TE_r(1))     , ') ms @Start'], ...
            ['Echo (', num2str(TE_r(ky_ctr)), ') ms @TE'   ], ...
            ['Echo (', num2str(TE_r(end))   , ') ms @End'  ], ...
            'Acquisition window' ) ;
end

hz = zoom( h0 ) ; 
set( hz, 'ActionPostCallback', { @cb_zoom_annot1, arg1 } ) ;











% ========================================================================
% Function to calculate FWHM
% - The FWHM is calcuated from the magnitude shape of the PSF, under the
%   assumption that the baseline signal (i.e. background noise) of the 
%   point-source is zero. 
% - The peak is normalized to 1 and the corresonding FWHM is returned.
% ========================================================================
function [x_fwhm, x_pos1, x_pos2] = get_fwhm_psf( ps_res ) 

[psf_max, xc_idx] = max( abs( ps_res ) ) ; 

vy = abs( ps_res )/psf_max ;
vx = 1:1:length( vy ) ;

xl_idx = find( vy( 1 : xc_idx-1 ) < 0.5 ) ;
xs1 = vx( xl_idx(end) ) ;
xs2 = vx( xl_idx(end)+1 ) ;
ys1 = vy( xl_idx(end) ) ;
ys2 = vy( xl_idx(end)+1 ) ;

xr_idx = find( vy( xc_idx+1 : end ) < 0.5 ) ;
xe1 = vx( xc_idx+xr_idx(1) ) ;
xe2 = vx( xc_idx+xr_idx(1)-1 ) ;
ye1 = vy( xc_idx+xr_idx(1) ) ;
ye2 = vy( xc_idx+xr_idx(1)-1 ) ;

x_pos1 = ( xs2 - xs1 )/( ys2 - ys1 )*( 0.5 - ys1 ) + xs1 ;
x_pos2 = ( xe2 - xe1 )/( ye2 - ye1 )*( 0.5 - ye1 ) + xe1 ;
x_fwhm  = x_pos2 - x_pos1 ;


% ========================================================================
% Function to calculate the normalized position
% ========================================================================
function [xn, yn] = get_norm_pos( h_axis, p1, p2 ) 

Xlm_h = h_axis.XLim ;
Ylm_h = h_axis.YLim ;
Pos_h = h_axis.Position ;
xn = Pos_h(1) + Pos_h(3)/( Xlm_h(2) - Xlm_h(1) )*( [p1(1) p2(1)] - Xlm_h(1) ) ;
yn = Pos_h(2) + Pos_h(4)/( Ylm_h(2) - Ylm_h(1) )*( [p1(2) p2(2)] - Ylm_h(1) ) ;


% ========================================================================
% Function to control figure objects
% ========================================================================
function cb_zoom_annot1( hObjet, eventdata, handles )

% Get current axis
hf_curr = gcf ;
ha_curr = hf_curr.CurrentAxes ;

% In case of the PSF plot
if ( strcmp( ha_curr.Tag, 'h2' ) ) ;

    % Delete
    harrow = findall( gcf, 'Tag', 'arrow_h2' ) ;
    for ii = 1 : 1 : size( harrow ) ;
        delete( harrow(ii) ) ;
    end

    % Get position
    pos_a   = handles{1}.Position ;
    pos_ax1 = pos_a(1) ; pos_ax2 = pos_a(1) + pos_a(3) ;
    pos_ay1 = pos_a(2) ; pos_ay2 = pos_a(2) + pos_a(4) ;

    % Arrow1
    [xn, yn] = get_norm_pos( handles{1}, handles{2}, handles{3} ) ;
    if ( ( yn(1) >= pos_ay1 ) && ( yn(1) <= pos_ay2 ) )

        if ( ( xn(1) >= pos_ax1 ) && ( xn(1) <= pos_ax2 ) && ( xn(2) >= pos_ax1 ) && ( xn(2) <= pos_ax2 ) )
            try
                ht = annotation( 'Arrow', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed1-1' ) ;
            end ;
        elseif ( ( xn(1) < pos_ax1 ) && ( xn(2) >= pos_ax1 ) && ( xn(2) <= pos_ax2 ) ) 
            xn(1) = pos_ax1 ; 
            try
                ht = annotation( 'Arrow', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed1-2' ) ;
            end ;
        elseif ( ( xn(1) >= pos_ax1 ) && ( xn(1) <= pos_ax2 ) && ( xn(2) > pos_ax2 ) ) 
            xn(2) = pos_ax2 ; 
            try
                ht = annotation( 'Line', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed1-3' ) ;
            end ;
        elseif ( ( xn(1) < pos_ax1 ) && ( xn(2) > pos_ax2 ) ) 
            xn(1) = pos_ax1 ; 
            xn(2) = pos_ax2 ; 
            try
                ht = annotation( 'Line', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed1-4' ) ;
            end ;
        end

    end 

    % Arrow2
    [xn, yn] = get_norm_pos( handles{1}, handles{4}, handles{5} ) ;
    if ( ( yn(1) >= pos_ay1 ) && ( yn(1) <= pos_ay2 ) )

        if ( ( xn(2) >= pos_ax1 ) && ( xn(2) <= pos_ax2 ) && ( xn(1) >= pos_ax1 ) && ( xn(1) <= pos_ax2 ) )
            try
                ht = annotation( 'Arrow', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed2-1' ) ;
            end ;
        elseif ( ( xn(2) >= pos_ax1 ) && ( xn(2) <= pos_ax2 ) && ( xn(1) > pos_ax2 ) ) 
            xn(1) = pos_ax2 ; 
            try
                ht = annotation( 'Arrow', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed2-2' ) ;
            end ;
        elseif ( ( xn(2) < pos_ax1 ) && ( xn(1) >= pos_ax1 ) && ( xn(1) <= pos_ax2 ) ) 
            xn(2) = pos_ax1 ; 
            try
                ht = annotation( 'Line', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed2-3' ) ;
            end ;
        elseif ( ( xn(2) < pos_ax1 ) && ( xn(1) > pos_ax2 ) ) 
            xn(2) = pos_ax1 ; 
            xn(1) = pos_ax2 ; 
            try
                ht = annotation( 'Line', xn, yn ) ;
                set( ht, 'Tag', 'arrow_h2' ) ;
            catch
                disp( 'Annotation drawing failed2-4' ) ;
            end ;
        end

    end 
    
% In case of the PSF plot
elseif ( strcmp( ha_curr.Tag, 'h4' ) ) ;

    % Delete
    harrow = findall( gcf, 'Tag', 'arrow_h4' ) ;
    for ii = 1 : 1 : size( harrow ) ;
        delete( harrow(ii) ) ;
    end

    % Get position
    pos_b   = handles{6}.Position ;
    pos_ax3 = pos_b(1) ; pos_ax4 = pos_b(1) + pos_b(3) ;

    % Arrow1
    p3 = handles{7} ; p3(2) = handles{6}.YLim(1) ;
    p4 = handles{8} ; p4(2) = handles{6}.YLim(1) ;
    [xm, ym] = get_norm_pos( handles{6}, p3, p4 ) ;
    if ( ( xm(1) >= pos_ax3 ) && ( xm(1) <= pos_ax4 ) && ( xm(2) >= pos_ax3 ) && ( xm(2) <= pos_ax4 ) )
        try
            ht = annotation( 'doublearrow', xm, ym, 'Color', [0 0 1] ) ;
            set( ht, 'Tag', 'arrow_h4' ) ;
        catch
            disp( 'Annotation drawing failed1-1' ) ;
        end ;
    elseif ( ( xm(1) < pos_ax3 ) && ( xm(2) >= pos_ax3 ) && ( xm(2) <= pos_ax4 ) ) 
        xm(1) = pos_ax3 ; 
        try
            ht = annotation( 'Arrow', xm, ym, 'Color', [0 0 1] ) ;
            set( ht, 'Tag', 'arrow_h4' ) ;
        catch
            disp( 'Annotation drawing failed1-2' ) ;
        end ;
    elseif ( ( xm(1) >= pos_ax3 ) && ( xm(1) <= pos_ax4 ) && ( xm(2) > pos_ax4 ) ) 
        xm(2) = pos_ax4 ; 
        try
            ht = annotation( 'Arrow', [xm(2) xm(1)], ym, 'Color', [0 0 1] ) ;
            set( ht, 'Tag', 'arrow_h4' ) ;
        catch
            disp( 'Annotation drawing failed1-3' ) ;
        end ;
    elseif ( ( xm(1) < pos_ax3 ) && ( xm(2) > pos_ax4 ) ) 
        xm(1) = pos_ax3 ; 
        xm(2) = pos_ax4 ; 
        try
            ht = annotation( 'Line', xm, ym, 'Color', [0 0 1] ) ;
            set( ht, 'Tag', 'arrow_h4' ) ;
        catch
            disp( 'Annotation drawing failed1-4' ) ;
        end ;
    end
    
    % Arrow2 (PF case)
    p3 = handles{9}  ; p3(2) = handles{6}.YLim(1) ;
    p4 = handles{10} ; p4(2) = handles{6}.YLim(1) ;
    if ( p4(1) > p3(1) )
        
        [xm, ym] = get_norm_pos( handles{6}, p3, p4 ) ;
        if ( ( xm(1) >= pos_ax3 ) && ( xm(1) <= pos_ax4 ) && ( xm(2) >= pos_ax3 ) && ( xm(2) <= pos_ax4 ) )
            try
                ht = annotation( 'doublearrow', xm, ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-1' ) ;
            end ;
        elseif ( ( xm(1) < pos_ax3 ) && ( xm(2) >= pos_ax3 ) && ( xm(2) <= pos_ax4 ) ) 
            xm(1) = pos_ax3 ; 
            try
                ht = annotation( 'Arrow', xm, ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-2' ) ;
            end ;
        elseif ( ( xm(1) >= pos_ax3 ) && ( xm(1) <= pos_ax4 ) && ( xm(2) > pos_ax4 ) ) 
            xm(2) = pos_ax4 ; 
            try
                ht = annotation( 'Arrow', [xm(2) xm(1)], ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-3' ) ;
            end ;
        elseif ( ( xm(1) < pos_ax3 ) && ( xm(2) > pos_ax4 ) ) 
            xm(1) = pos_ax3 ; 
            xm(2) = pos_ax4 ; 
            try
                ht = annotation( 'Line', xm, ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-4' ) ;
            end ;
        end
        
    end 
    
    % Arrow2 (PF 2-side case)
    p3 = handles{11} ; p3(2) = handles{6}.YLim(1) ;
    p4 = handles{12} ; p4(2) = handles{6}.YLim(1) ;
    if ( p4(1) > p3(1) )
        
        [xm, ym] = get_norm_pos( handles{6}, p3, p4 ) ;
        if ( ( xm(1) >= pos_ax3 ) && ( xm(1) <= pos_ax4 ) && ( xm(2) >= pos_ax3 ) && ( xm(2) <= pos_ax4 ) )
            try
                ht = annotation( 'doublearrow', xm, ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-1' ) ;
            end ;
        elseif ( ( xm(1) < pos_ax3 ) && ( xm(2) >= pos_ax3 ) && ( xm(2) <= pos_ax4 ) ) 
            xm(1) = pos_ax3 ; 
            try
                ht = annotation( 'Arrow', xm, ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-2' ) ;
            end ;
        elseif ( ( xm(1) >= pos_ax3 ) && ( xm(1) <= pos_ax4 ) && ( xm(2) > pos_ax4 ) ) 
            xm(2) = pos_ax4 ; 
            try
                ht = annotation( 'Arrow', [xm(2) xm(1)], ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-3' ) ;
            end ;
        elseif ( ( xm(1) < pos_ax3 ) && ( xm(2) > pos_ax4 ) ) 
            xm(1) = pos_ax3 ; 
            xm(2) = pos_ax4 ; 
            try
                ht = annotation( 'Line', xm, ym, 'Color', [1 0 0] ) ;
                set( ht, 'Tag', 'arrow_h4' ) ;
            catch
                disp( 'Annotation drawing failed1-4' ) ;
            end ;
        end
        
    end 

else
end ;

